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
                 'HA': ha,
                 'h_acc': h_acceptors,
                 'h_don': h_donors,
                 'rings': rings
                 }

    return desc_dict

def make_r_sprites(r_group, label):
    """Function to get all the unique smiles for a given r group tag and make the .pngs. Assumes the smiles are sorted
    i.e. that the tags e.g. A01 increases monotonically.
    """
    data = pd.read_csv('../data/r_group_decomp.csv')
    r_smiles = data[r_group].unique()
    for i, r in enumerate(r_smiles):
        tag = data[label].unique()[i]
        mol = Chem.MolFromSmiles(r)
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(mol)
        d.FinishDrawing()
        d.WriteDrawingText(f'../Images/r_group_pngs/{tag}.png')

def lipinski(desc_dict):
    """Calculate Lipinski from the descriptor dictionary.
    Return the number of rules broken and whether the molecule passes.
    """
    violations = 0
    if desc_dict['MW'] > 500: violations += 1
    if desc_dict['h_acc'] > 10: violations += 1
    if desc_dict['h_don'] > 5: violations += 1
    if desc_dict['logP'] > 5: violations += 1
    if violations > 1:
        result = 'fails'
    else:
        result = 'passes'
    
    return violations, result

# make_r_sprites('R2','btag')
