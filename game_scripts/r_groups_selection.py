from rdkit import Chem
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.Chem import Draw, AllChem
import pandas as pd


def get_selection(r_group, no_picks, data):
    """Selects a specified number of either R1 or R2 groups from a provided Dataframe of groups.

    :param r_group: either 'R1' or 'R2'
    :type r_group: string
    :param no_picks: number of R groups to retrieve
    :type no_picks: int
    :param data: Pandas dataframe containing all the R1 and R2 groups
    :type data: :class:`pandas.Dataframe`
    :return: [[R group identifiers], [R group SMILES]]
    :rtype: nested list
    """
    r_groups = data[r_group.upper()]  # subset the data to only the given r group (either 'R1' or 'R2')
    mols = [Chem.MolFromSmiles(mol) for mol in r_groups]
    fps = [Chem.AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in mols]  # make fingerprints
    picker = rdSimDivPickers.MaxMinPicker()
    picks = list(picker.LazyBitVectorPick(fps, len(fps), no_picks))
    if r_group == 'R1':
        return list(data['atag'][picks]), list(data['R1'][picks])
    elif r_group == 'R2':
        return list(data['btag'][picks]), list(data['R2'][picks])
