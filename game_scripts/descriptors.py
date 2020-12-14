from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Crippen
import pandas as pd


def get_descriptors(mol):
    """Function to get the specified RDkit descriptors for a given molecule. Gets the Lipinski parameters.
    Takes a SMILE string, converts to molecule and then calculates.

    :param mol: the molecule to calculate descriptors for
    :type mol: SMILE string

    :return: dictionary with descriptor name-value pairs for mol (smile), molecular weight, logP, TPSA,
    Heavy Atom count, h_acceptors, h_donors, ring count.
    :rtype: dict
    """
    smile = mol
    mol = Chem.MolFromSmiles(mol)

    # Calculate descriptors
    mw = Descriptors.ExactMolWt(mol)
    log_p = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)  # topological polar surface area
    ha = Lipinski.HeavyAtomCount(mol)  # heavy atom count
    h_acceptors = Lipinski.NumHAcceptors(mol)
    h_donors = Lipinski.NumHDonors(mol)
    rings = Lipinski.RingCount(mol)

    desc_dict = {'mol': smile,
                 'MW': f'{mw:.4f}',
                 'logP': log_p,
                 'TPSA': tpsa,
                 'Heavy Atoms': ha,
                 'h_acceptors': h_acceptors,
                 'h_donors': h_donors,
                 'rings': rings
                 }

    return desc_dict


x = get_descriptors(mol='Oc1ccc(C[*:1])cc1')
print(x)