import copy
import socket

import parsl
from parsl.config import Config
from parsl.executors import ThreadPoolExecutor

try:
    from SocketServer import ThreadingMixIn
except ImportError:
    from socketserver import ThreadingMixIn

import numpy as np

from openeye import oechem, oedocking

from xmlrpc.client import Binary

from xmlrpc.server import SimpleXMLRPCServer, SimpleXMLRPCRequestHandler

import argparse

class OEDockingServer:
    def __init__(self, receptors, names):
        self.receptors = {}
        self.results = {}
        self.idx = 0

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

        self.idx += 1
        if isinstance(smiles, list):
            self.results[self.idx] = []
            for smile in smiles :
                idx_ =  oedock_from_smiles(receptor, smile, oe_options=oe_options)
                self.results[self.idx].append(idx_)
                print("sent", smile, idx_, self.idx)
        else:
            idx_ =  oedock_from_smiles(receptor, smiles, oe_options=oe_options)
            self.results[self.idx] = idx_
            print("sent", smiles, idx_, self.idx)

        return self.idx

    def QueryStatus(self, queryidx):
        if not isinstance(self.results[queryidx], list):
            print("no list")
            return self.results[queryidx].done()
        else:
            done = True
            for res in self.results[queryidx]:
                done = done and res.done()
            return done

    def QueryResults(self, queryidx):
        if not isinstance(self.results[queryidx], list):
            print("no list")
            results = self.results[queryidx].result()
        else:
            results = []
            for res in self.results[queryidx]:
                results.append(res.result())
            print(results)
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
