from rdkit import Chem
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.Chem import Draw, AllChem
import pandas as pd


def get_selection(r_group, no_picks, data):
    r_groups = data[r_group]
    mols = [Chem.MolFromSmiles(mol) for mol in r_groups]
    fps = [Chem.AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in mols]
    picker = rdSimDivPickers.MaxMinPicker()
    picks = list(picker.LazyBitVectorPick(fps, len(fps), no_picks))
    mol_picks = [mols[ind] for ind in picks]
    if r_group == 'R1':
        return list(data['atag'][picks]), list(data['R1'][picks])
    elif r_group == 'R2':
        return list(data['btag'][picks]), list(data['R2'][picks])
