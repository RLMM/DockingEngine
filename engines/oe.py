import numpy as np
from openeye import oechem, oeomega, oedocking, oequacpac
from parsl import python_app


# Molecule Generation
def mol_from_smiles(smiles):
    mol = oechem.OEMol()
    if not oechem.OESmilesToMol(mol, smiles):
        print("SMILES invalid for string", smiles)
        raise ValueError("SMILES invalid for string", smiles)
    return mol

class OEOptions:
    def __init__(self, **kwargs):
        self.num_sterocenters = 6
        self.force_flipper = False
        self.use_flipper = True
        self.use_tautomer = True

        self.use_hybrid = True
        self.high_resolution = True
        self.cache_receptor = True
        self.cache_receptor_clear = False

        self.num_poses = 1
        self.nan_to_none = True

        self.update(kwargs)

    def update(self, d : dict):
        for k, v in d.items():
            if k in self.__dict__:
                self.__dict__[k] = v
            else:
                print("not in oeoptions", k)

# Conformer Generation
def enumerate_from_smiles(smiles, oe_options=None):
    oe_options = oe_options or OEOptions()

    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Pose)
    omegaOpts.SetMaxSearchTime(60)
    omega = oeomega.OEOmega(omegaOpts)

    if oe_options.use_tautomer:
        tautomer_options = oequacpac.OETautomerOptions()
        pKa_norm = True
        taut_iter = lambda x: oequacpac.OEGetReasonableTautomers(x, tautomer_options, pKa_norm)
    else:
        taut_iter = lambda x: [x]

    if oe_options.use_flipper:
        flipper = lambda x: oeomega.OEFlipper(x.GetActive(), oe_options.num_sterocenters, oe_options.force_flipper)
    else:
        flipper = lambda x: [x]

    # get molecule
    try:
        molecule = mol_from_smiles(smiles)
    except ValueError:
        return [None]

    results = []
    for enantiomer in flipper(molecule):
        for tautomer in taut_iter(enantiomer):
            tautomer = oechem.OEMol(tautomer)
            omega.Build(tautomer)
            tautomer2 = oechem.OEMol(tautomer)
            results.append(tautomer2)
    return results

def init_oedock_from_receptor(receptor, oe_options=None):
    oe_options = oe_options or OEOptions()

    dockMethod = oedocking.OEDockMethod_Hybrid if oe_options.use_hybrid else oedocking.OEDockMethod_Chemgauss4
    dockResolution = oedocking.OESearchResolution_High if oe_options.high_resolution else oedocking.OESearchResolution_Default
    dock = oedocking.OEDock(dockMethod, dockResolution)
    assert (dock.Initialize(receptor))

    if (not oedocking.OEReceptorHasCachedScore(receptor)) and oe_options.cache_receptor:
        dock.CacheScoringSetup(receptor, clearOldData=oe_options.cache_receptor_clear)
    return dock

def receptor_from_file(receptor_filename):
    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, receptor_filename):
        oechem.OEThrow.Fatal("Unable to read receptor")

    if not oedocking.OEReceptorHasBoundLigand(receptor):
        raise Exception("Receptor does not have bound ligand")
    return receptor

def setup_receptor_from_file(receptor_filename, oe_options=None):
    oe_options = oe_options or OEOptions()

    receptor = receptor_from_file(receptor_filename)
    dock = init_oedock_from_receptor(receptor, oe_options)

    return dock, receptor


def dock_molecule(smiles, dock_obj, oe_options=None):
    oe_options = oe_options or OEOptions()

    scores = []
    confs = enumerate_from_smiles(smiles, oe_options=oe_options)
    for conf in confs:
        if conf is None:
            continue
        score = safe_dock_score_(dock_obj, conf)
        if score is not None:
            scores.append(score)

    if len(scores) > 0:
        return min(scores)
    else:
        return np.nan


def safe_dock_score_(dock, mol, oe_options=None):
    oe_options = oe_options or OEOptions()

    try:
        dockedMol = oechem.OEMol()
        retCode = dock.DockMultiConformerMolecule(dockedMol, mol, oe_options.num_poses)
        if (retCode != oedocking.OEDockingReturnCode_Success):
            rets = oedocking.OEDockingReturnCodeGetName(retCode)
            print('oeodock error', rets)
            return None
        score = dock.ScoreLigand(dockedMol)
        if oe_options.nan_to_none and (score > 10000 or np.isnan(score) or np.isinf(score)):
            return None
        else:
            return score
    except:
        return None


@python_app
def oedock_from_smiles(receptor, smiles, oe_options=None):
    oe_options = oe_options or OEOptions()

    docker = init_oedock_from_receptor(receptor, oe_options=oe_options)

    if isinstance(smiles, list):
        scores = []
        for smile in smiles:
            scores.append(dock_molecule(smile, docker, oe_options=oe_options))
        return scores
    else:
        score = dock_molecule(smiles, docker, oe_options=oe_options)
        return score
