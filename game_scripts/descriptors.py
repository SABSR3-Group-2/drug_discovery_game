from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
from rdkit.Chem import Crippen
import pandas as pd
"""Feedback from Roche on 16/12/2020: sort R groups by the basic features calculated in this script. When/if
sorting the final compounds you'd want to assay them by more complex traits such as lipopilicity, 
metabolism etc. (the ones given in the presentation they shared"""

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

def make_r_sprites(r_group, label):
    """Function to get all the unique smiles for a given r group tag and make the .pngs.

    """
    data = pd.read_csv('../data/r_group_decomp.csv')
    r_smiles = data[r_group].unique()
    for i, r in enumerate(r_smiles):
        mol = Chem.MolFromSmiles(r)
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(mol)
        d.FinishDrawing()
        d.WriteDrawingText(f'../Images/{label}{i + 1}.png')
