from __future__ import print_function

import argparse
import sys
import time
from xmlrpc.client import ServerProxy, Binary, Fault

from openeye import oechem, oedocking


class DockingClient:
    def __init__(self, hostname, port):
        self.hostname = hostname
        self.port = port

    def __call__(self, smiles, receptor_name, receptor=None, receptor_path=None):
        receptor = self.get_receptor_from_request(receptor, receptor_path)
        s = ServerProxy("http://" + f'{self.hostname}:{self.port}')
        idx = s.SubmitQuery(smiles, receptor, receptor_name)

        not_done = True
        while not_done:
            try:
                not_done = not s.QueryStatus(idx)
                time.sleep(2)
            except Fault as e:
                print(str(e), file=sys.stderr)
                return 1

        results = s.QueryResults(idx)
        return results

    def get_receptor_from_request(self, receptor=None, receptor_path=None):
        if receptor is not None and receptor_path is not None:
            print("can only use one")
            assert (False)

        if receptor_path is not None:
            return receptor_path
        else:
            bytes = oedocking.OEWriteReceptorToBytes('.oeb', receptor)
            return Binary(bytes)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', help='port for server', type=int, default=8080)
    parser.add_argument('-a', help='server portname', type=str, default='0.0.0.0')
    parser.add_argument('-s', help='smiles', type=str, default='CC(C)[C@H](C(=O)[O-])NC(=O)CCCCSc1c2c([nH]cn2)nc(n1)N')
    parser.add_argument('-r', help='receptor', type=str, default='/Users/austin/gly_adrp.oeb')
    parser.add_argument('-n', help='receptor_name', type=str, default=None)
    return parser.parse_args()


def main(args):
    smiles = args.s
    receptor_name = args.r
    hostname = args.a
    port = args.p

    receptor = oechem.OEGraphMol()
    oedocking.OEReadReceptorFile(receptor, args.r)
    client = DockingClient(hostname, port)
    client(smiles, receptor_name, receptor=receptor)

    return 0


if __name__ == '__main__':
    args = get_args()
    main(args)
