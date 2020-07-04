import socket

from parsl.config import Config
from parsl.executors import ThreadPoolExecutor

try:
    from SocketServer import ThreadingMixIn
except ImportError:
    from socketserver import ThreadingMixIn

from openeye import oechem, oedocking

from xmlrpc.client import Binary
from xmlrpc.server import SimpleXMLRPCServer, SimpleXMLRPCRequestHandler

from engines.oe import oedock_from_smiles, setup_receptor_from_file, OEOptions
import argparse
import parsl

class OEDockingServer:
    def __init__(self, receptors, names):
        self.receptors = {}
        self.results = {}
        self.idx = 0

        assert (len(receptors) == len(names))
        for res in zip(receptors, names):
            self.get_receptor(*res)

    def get_receptor(self, receptor, name):
        self.receptors = {}
        if isinstance(receptor, Binary):
            receptor_mol = oechem.OEGraphMol()
            oedocking.OEReadReceptorFromBytes(receptor_mol, '.oeb', receptor.data)
            self.receptors[name] = receptor_mol
        elif name not in self.receptors:
            _, receptor = setup_receptor_from_file(receptor)
            self.receptors[name] = receptor
        return self.receptors[name]

    def SubmitQuery(self, smiles, receptor, receptor_name, oe_options=None):
        if oe_options is None:
            oe_options = OEOptions()
        else:
            oe_options = OEOptions(**oe_options)

        receptor = self.get_receptor(receptor, receptor_name)
        self.idx += 1
        self.results[self.idx] = oedock_from_smiles(receptor, smiles, oe_options=oe_options)
        return self.idx

    def QueryStatus(self, queryidx):
        return self.results[queryidx].done()

    def QueryResults(self, queryidx):
        results = self.results[queryidx].result()
        del self.results[queryidx]
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

    parsl.load()
    print("Parsl loaded.")
    sender(args.a, args.p, args.receptors, args.names)
