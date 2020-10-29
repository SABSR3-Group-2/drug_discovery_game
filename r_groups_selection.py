from rdkit import Chem
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.Chem import Draw, AllChem
import pandas as pd
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

def get_selection(r_group, no_picks, DataSource='data/r_groups_pIC50.csv):
    data = pd.read_csv(DataSource, delimiter=',')
    r_groups = data[r_group]
    mols = [Chem.MolFromSmiles(mol) for mol in r_groups]
    fps = [Chem.AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in mols]
    picker = rdSimDivPickers.MaxMinPicker()
    picks = list(picker.LazyBitVectorPick(fps, len(fps), no_picks))
    mol_picks = [mols[ind] for ind in picks]
    
    if r_group == 'R1':
        return Draw.MolsToGridImage(mol_picks, legends=list(data['Atag'][picks]))

    elif r_group == 'R2':
        return Draw.MolsToGridImage(mol_picks, legends=list(data['Btag'][picks]))
