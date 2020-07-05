import copy
import socket
import time
from concurrent.futures import wait, as_completed, FIRST_COMPLETED
import parsl
from parsl.config import Config
from parsl.executors import ThreadPoolExecutor

try:
    from SocketServer import ThreadingMixIn
except ImportError:
    from socketserver import ThreadingMixIn
from threading import Thread
from threading import Lock
from threading import Condition
from threading import Event
from openeye import oechem, oedocking

from xmlrpc.client import Binary

from xmlrpc.server import SimpleXMLRPCServer, SimpleXMLRPCRequestHandler

import argparse


class ReadWriteLock(object):
    """ Basic locking primitive that allows multiple readers but only
    a single writer at a time. Useful for synchronizing database
    updates. Priority is given to pending writers. """
    def __init__(self):
        self.cond = Condition()
        self.readers = 0
        self.writers = 0

    def AcquireReadLock(self):
        self.cond.acquire()
        try:
            while self.writers:
                self.cond.wait()

            self.readers += 1
            assert self.writers == 0
        finally:
            self.cond.notifyAll()
            self.cond.release()

    def ReleaseReadLock(self):
        self.cond.acquire()
        assert self.readers > 0
        try:
            self.readers -= 1
        finally:
            self.cond.notifyAll()
            self.cond.release()

    def AcquireWriteLock(self):
        self.cond.acquire()
        self.writers += 1
        while self.readers:
            self.cond.wait()

        assert self.readers == 0
        assert self.writers > 0

    def ReleaseWriteLock(self):
        assert self.readers == 0
        assert self.writers > 0

        self.writers -= 1
        self.cond.notifyAll()
        self.cond.release()

class OEDockingServer:
    def __init__(self, receptors, names):
        self.receptors = {}
        self.results = {}
        self.idx = 0
        self.done_arr = {}
        self.waiting_futures = {}
        self.query_time = {}
        self.lock = Lock()

        assert (len(receptors) == len(names))
        for res in zip(receptors, names):
            self.get_receptor(*res)

    def AddReceptor(self, receptor, name):
        if isinstance(receptor, Binary):
            receptor_mol = oechem.OEGraphMol()
            oedocking.OEReadReceptorFromBytes(receptor_mol, '.oeb', receptor.data)
            self.receptors[name] = oechem.OEGraphMol(receptor_mol)
        else:
            _, receptor_mol = setup_receptor_from_file(receptor)
            self.receptors[name] = oechem.OEGraphMol(receptor_mol)

    def get_receptor(self, receptor, name):
        if name in self.receptors:
            return oechem.OEGraphMol(self.receptors[name])
        elif isinstance(receptor, Binary):
            receptor_mol = oechem.OEGraphMol()
            oedocking.OEReadReceptorFromBytes(receptor_mol, '.oeb', receptor.data)
            self.receptors[name] = receptor_mol
            return receptor_mol
        else:
            _, receptor_mol = setup_receptor_from_file(receptor)
            self.receptors[name] = receptor_mol
            return receptor_mol

    def SubmitQuery(self, smiles, receptormol, receptor_name, oe_options=None, chunksize=1):
        receptor = self.get_receptor(receptormol, receptor_name)
        if oe_options is None:
            oe_options = OEOptions()
        else:
            oe_options = OEOptions(**oe_options)

        self.lock.acquire()
        try:
            self.idx += 1
            cur_idx = self.idx
        except:
            exit()
        finally:
            self.lock.release()

        if isinstance(smiles, list):
            results_ = []
            done_arr_ = []
            waiting_futures_ = []

            for i, smile in enumerate(smiles):
                idx_ =  oedock_from_smiles(receptor, smile, oe_options=oe_options)
                results_.append(idx_)
                done_arr_.append(i)
                waiting_futures_.append(False)

            self.results[cur_idx] = results_
            self.done_arr[cur_idx] = done_arr_
            self.waiting_futures[cur_idx] = waiting_futures_
            print(f"[{time.time()}] queued {len(smiles)} smiles with job id {self.idx}")
        else:
            idx_ =  oedock_from_smiles(receptor, smiles, oe_options=oe_options)
            self.results[cur_idx] = idx_
            print(f"[{time.time()}] queued {smiles} with job id {self.idx}")

        self.query_time[cur_idx] = time.time()
        return self.idx

    def wait_for_change(self, queryidx):
        if not isinstance(self.results[queryidx], list):
            wait([self.results[queryidx]])
        else:
            futures = self.results[queryidx]
            futures_done, futures_not_done = wait([futures[i] for i in self.waiting_futures[queryidx]],
                                                  return_when=FIRST_COMPLETED)
            new_waitlist = []
            for j in range(len(futures)):
                if futures[j] in futures_not_done:
                    new_waitlist.append(j)
            self.waiting_futures[queryidx] = new_waitlist

    # def QueryStatus(self, queryidx, blocking=True):
    #     if blocking:
    #         self.wait_for_change(queryidx)
    #
    #     dcount = 0
    #     if not isinstance(self.results[queryidx], list):
    #         dcount = int(self.results[queryidx].done())
    #         total = 1
    #     else:
    #         total = len(self.results[queryidx])
    #         for i, res in enumerate(self.results[queryidx]):
    #             if self.done_arr[queryidx][i]:
    #                 is_done_ = True
    #             else:
    #                 is_done_ = res.done()
    #                 self.done_arr[queryidx][i] = is_done_
    #             dcount += int(is_done_)
    #
    #     if dcount == total:
    #         self.query_time[queryidx] = time.time() - self.query_time[queryidx]
    #
    #     return dcount, total

    def QueryStatus(self, queryidx, blocking=False):
        # if blocking:
        #     self.wait_for_change(queryidx)

        query = self.results[queryidx]
        if not isinstance(query, list):
            if not query.done():
                return False
        else:
            for i, res in enumerate(query):
                if not res.done():
                    return False

        qtime = self.query_time[queryidx]
        self.query_time[queryidx] = time.time() - qtime
        print(f"[{time.time()}] finished results for job {queryidx}.")

        return True

    def QueryResults(self, queryidx):
        results_queryidx = self.results[queryidx]
        if not isinstance(results_queryidx, list):
            results = results_queryidx.result()
            lenres = 1
        else:
            results = []
            for res in results_queryidx:
                results.append(res.result())
            lenres = len(results)
        print(f"[{time.time()}] sent results of size {lenres} for job {queryidx}. Total time {self.query_time[queryidx]}, {self.query_time[queryidx] / lenres}")

        return results


class RequestHandler(SimpleXMLRPCRequestHandler):
    rpc_paths = ('/RPC2',)


class AsyncXMLRPCServer(ThreadingMixIn, SimpleXMLRPCServer):
    # if a shutdown request occurs through a signal force everything to terminate immediately
    daemon_threads = True
    allow_reuse_address = True


def sender(hostname, portnumber, receptors, names):
    # create server
    server = AsyncXMLRPCServer((hostname, portnumber),
                               requestHandler=RequestHandler,
                               logRequests=False)

    hs_name = hostname
    hostname, portnumber = server.socket.getsockname()

    if hostname == "0.0.0.0":
        hostname = socket.gethostname()
        hs_name = 'localhost'

    print("Listening for DockingServer.py requests on %s (%s):%i\n\n" % (hostname, hs_name, portnumber))

    # register the XMLRPC methods
    server.register_introspection_functions()
    server.register_instance(OEDockingServer(receptors, names))

    try:
        server.serve_forever()
    finally:
        server.server_close()

    return 0


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', help='hostname', required=False, type=str, default='localhost')
    parser.add_argument('-p', help='port', required=False, type=int, default=8080)
    parser.add_argument('--receptors', action='append', type=str, required=False, help='receptor files to preload must just file path. Can provide many times')
    parser.add_argument('--named_receptors', action='append', type=str, required=False, help='receptor files to preload must be filepath:name. Can provide many times')
    parser.add_argument('--n_jobs', type=int, default=8, help='number of executors', required=False)
    args = parser.parse_args()

    receptors = []
    names = []

    if args.receptors is not None:
        for receptor in args.receptors:
            names.append(receptor.split('/')[-1].split('.')[0])
            receptors.append(receptor)
            print(f"Adding receptor {receptors[-1]}, and naming it {names[-1]}")

    if args.named_receptors is not None:
        for named_receptor in args.named_receptors:
            receptor, name = named_receptor.split(":")
            receptors.append(receptor)
            names.append(name)
            print(f"Adding receptor {receptors[-1]}, named {names[-1]}")

    args.receptors = receptors
    args.names = names
    return args


if __name__ == '__main__':

    args = get_args()

    config = Config(
        executors=[
            ThreadPoolExecutor(max_threads=args.n_jobs)
        ],
    )

    print("Parsl loaded.")

    from engines.oe import oedock_from_smiles, setup_receptor_from_file, OEOptions

    parsl.load(config)

    sender(args.a, args.p, args.receptors, args.names)
