import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, GraphDescriptors

test_molecule = ['mjgklghsdjk','CCCCCC','C(Cl)Cl']


def desc_calculation(molecules):
    """Calculates 119 RDKit descriptors for one molecule or list of molecules"""
    mtype = type(molecules)

    if mtype == list:
        return _set_desc_calculeation(molecules)
    elif mtype == str:
        return _set_desc_calculeation([molecules])
    return None


def _set_desc_calculeation(molecules):
    set_descriptors = list()
    error_mole = list()
    error_mole_index = list()
    for i in range(len(molecules)):
        try:
            set_descriptors.append(_compound_desc_calculation(molecules[i]))
        except Exception as e:
            if len(molecules) == 1:
                return None, e
            else:
                error_mole.append(molecules[i])
                error_mole_index.append(i)

    return np.asarray(set_descriptors), error_mole_index, error_mole


def _compound_desc_calculation(smiles):
    list_desc = ['MolLogP', 'MolMR', 'LabuteASA', 'TPSA', 'MolWt', 'ExactMolWt', 'CalcNumLipinskiHBA',
                 'CalcNumLipinskiHBD', 'NumRotatableBonds',
                 'CalcNumHBD', 'CalcNumHBA', 'CalcNumAmideBonds', 'NumHeteroatoms', 'HeavyAtomCount', 'NumAtoms',
                 'CalcNumAtomStereoCenters',
                 'CalcNumUnspecifiedAtomStereoCenters', 'RingCount', 'NumAromaticRings', 'NumSaturatedRings',
                 'NumAliphaticRings',
                 'NumAromaticHeterocycles', 'NumSaturatedHeterocycles', 'NumAliphaticHeterocycles',
                 'NumAromaticCarbocycles',
                 'NumSaturatedCarbocycles', 'NumAliphaticCarbocycles', 'FractionCSP3',
                 'Chi0v', 'Chi1v', 'Chi2v', 'Chi3v', 'Chi4v', 'Chi1n', 'Chi2n', 'Chi3n', 'Chi4n',
                 'HallKierAlpha', 'Kappa1', 'Kappa2', 'Kappa3',
                 'SlogP_VSA1', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7',
                 'SlogP_VSA8', 'SlogP_VSA9', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12',
                 'SMR_VSA1', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8',
                 'SMR_VSA9', 'SMR_VSA10',
                 'PEOE_VSA1', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8',
                 'PEOE_VSA9', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14',
                 'MQNs_']
    calculated_descriptors = list()
    rdkit_mol = Chem.MolFromSmiles(smiles)

    if rdkit_mol is None:
        raise Exception('bad molecule {}!'.format(smiles))

    Chem.Kekulize(rdkit_mol)
    for i in range(len(list_desc)):
        feature = list_desc[i]
        if feature == 'NumAtoms':
            tmp = Chem.rdmolops.AddHs(rdkit_mol)
            look = len(Chem.rdchem.Mol.GetAtoms(tmp))
        else:
            # f = lambda descrs: [getattr(field, feature)(rdkit_mol) for field in descrs if hasattr(field, feature)]
            # look = f([Descriptors, rdMolDescriptors, GraphDescriptors])[0]
            descrs = [Descriptors, rdMolDescriptors, GraphDescriptors]
            look = next((getattr(field, feature)(rdkit_mol) for field in descrs if hasattr(field, feature)), None)
        if type(look) == list:
            calculated_descriptors.extend(look)
        else:
            calculated_descriptors.append(look)

    calculated_descriptors = np.asarray(calculated_descriptors)
    return calculated_descriptors
