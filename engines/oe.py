import numpy as np
from openeye import oechem, oeomega, oedocking, oequacpac
from parsl import python_app


# Molecule Generation
def mol_from_smiles(smiles):
    mol = oechem.OEMol()
    if not oechem.OESmilesToMol(mol, smiles):
        raise ValueError("SMILES invalid for string", smiles)
    else:
        return mol


# Conformer Generation
def enumerate_from_smiles(smiles, num_sterocenters=6, force_flipper=False, use_flipper=True, use_tautomer=False):
    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Pose)
    omega = oeomega.OEOmega(omegaOpts)

    if use_tautomer:
        tautomer_options = oequacpac.OETautomerOptions()
        pKa_norm = True
        taut_iter = lambda x: oequacpac.OEGetReasonableTautomers(x, tautomer_options, pKa_norm)
    else:
        taut_iter = lambda x: [x]

    if use_flipper:
        flipper = lambda x: oeomega.OEFlipper(x.GetActive(), num_sterocenters, force_flipper)
    else:
        flipper = lambda x: [x]

    # get molecule
    try:
        molecule = mol_from_smiles(smiles)
    except ValueError:
        print("NAN")
        return [None]

    results = []
    for enantiomer in flipper(molecule):
        for tautomer in taut_iter(enantiomer):
            tautomer = oechem.OEMol(tautomer)
            omega.Build(tautomer)
            tautomer2 = oechem.OEMol(tautomer)
            results.append(tautomer2)
    return results


# Docking Engine
def from_cache(receptor):
    print("in from cache")
    dockMethod = oedocking.OEDockMethod_Hybrid
    dockResolution = oedocking.OESearchResolution_High
    dock = oedocking.OEDock(dockMethod, dockResolution)
    if not dock.Initialize(receptor):
        print("NO WQORKING!!")
    print("out from cache")
    return dock


def get_dock_obj(receptor_filename):
    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, receptor_filename):
        oechem.OEThrow.Fatal("Unable to read receptor")

    if not oedocking.OEReceptorHasBoundLigand(receptor):
        raise Exception("Receptor does not have bound ligand")

    dockMethod = oedocking.OEDockMethod_Hybrid
    dockResolution = oedocking.OESearchResolution_High
    dock = oedocking.OEDock(dockMethod, dockResolution)
    assert (dock.Initialize(receptor))
    dock.CacheScoringSetup(receptor, clearOldData=True)
    return dock, receptor


def dock_molecule(smiles, dock_obj, num_sterocenters=12, force_flipper=False, use_flipper=False, use_tautomer=True):
    scores = []
    print(smiles, dock_obj)
    confs = enumerate_from_smiles(smiles, use_flipper=use_flipper, use_tautomer=use_tautomer,
                                  num_sterocenters=num_sterocenters, force_flipper=force_flipper)
    for conf in confs:
        print(conf)
        if conf is None:
            continue
        score = safe_dock_score_(dock_obj, conf)
        if score is not None:
            scores.append(score)

    if len(scores) > 0:
        return min(scores)
    else:
        return 123092.0


def safe_dock_score_(dock, mol, num_poses=1, nan_to_none=True):
    try:
        dockedMol = oechem.OEMol()
        retCode = dock.DockMultiConformerMolecule(dockedMol, mol, num_poses)
        if (retCode != oedocking.OEDockingReturnCode_Success):
            rets = oedocking.OEDockingReturnCodeGetName(retCode)
            print(rets)
            return None
        score = dock.ScoreLigand(dockedMol)
        if nan_to_none and (score > 10000 or np.isnan(score) or np.isinf(score)):
            return None
        else:
            return score
    except:
        return None


@python_app
def oedock_from_smiles(receptor, smiles, num_stero=4, use_tautomers=True, force_flipper=False, use_flipper=True):
    docker = from_cache(receptor)
    score = dock_molecule(smiles, docker, num_sterocenters=4, use_tautomer=True, force_flipper=False, use_flipper=True)
    print(score)
    return score
