from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem import AllChem
from rdkit import RDLogger
from rdkit.SimDivFilters import rdSimDivPickers
import pandas as pd


def remove_h(r):
    """ Removes unwanted hydrogens and sanitizes mol

    :param r: an r group for H removal
    :type r: smile

    :return: sanitized r group as a Smile
    :rtype: Smile string
    """
    RDLogger.DisableLog('rdApp.*')  # don't print warning message for hydrogen removal
    if r == '[H][*:1].[H][*:1]':
        return '[H][*:1]'
    else:
        mol = Chem.MolFromSmiles(r)
        substruct = Chem.MolFromSmiles('[H][*:1]')
        mol = AllChem.DeleteSubstructs(mol, substruct) # remove hydrogen
        Chem.SanitizeMol(mol) # sanitize resulting molecule
        return Chem.MolToSmiles(mol)


def read_mols(filename, scaffold):
    """Reads the original csv and extracts the molecules, tags and pIC50 values. With large files the
    R group decomposition step can take a long time. Recommended to write the output to a file.

    :param filename: path to the .csv file containing the data e.g. `data/pickett_data.csv`
    :type filename: str
    :param scaffold: the core of the molecules as a Smile
    :type scaffold: Smile string

    :return: :class:`pandas.DataFrame` with columns ['mol','core','atag','R1','btag','R2',etc...,'pic50']
    :rtype: :class:`pandas.DataFrame`
    """
    data = pd.read_csv(filename, delimiter=',')
    data.columns = map(str.lower, data.columns)
    scaffold = Chem.MolFromSmiles(scaffold)
    mols = [Chem.MolFromSmiles(x) for x in data['smiles']]  # get all the mols from the smiles
    groups, _ = rdRGD.RGroupDecompose([scaffold], mols, asSmiles=True) # run r group decomp
    groups_frame = pd.DataFrame(groups)
    R_cols = [col for col in groups_frame.columns if 'R' in col]
    for col in R_cols: # add tag identifier to r groups e.g. A01
        num = int(col[1])
        letter = chr(ord('a') + num - 1)
        tag = letter + 'tag'
        groups_frame.insert(
            loc=groups_frame.columns.get_loc(col),
            column=tag,
            value=groups_frame[col].factorize()[0] + 1
        )
        groups_frame[tag] = groups_frame[tag].apply("{:02d}".format)
        groups_frame[tag] = letter.upper() + groups_frame[tag].astype(str)
        groups_frame[col] = groups_frame[col].apply(remove_h) # remove unwanted hydrogens
    # add columns for assay data and molecules
    groups_frame['pic50'] = data['pic50_mmp12']
    groups_frame['clearance_mouse'] = data['mlclearance (hep,m3c)']
    groups_frame['clearance_human'] = data['mlclearance (mic,h3c)']
    groups_frame['logd'] = data['mllogd']
    groups_frame['pampa'] = data['mlpampa (2c)']
    groups_frame.insert(0, 'mol', mols)

    return groups_frame


def rank_mols(data, feature):
    """Ranks (or clusters) the molecules according to specified feature which must a column name of data
    is pIC50

    :param data: dataframe with processed data. Outputted by read_mols
    :type data: :class:`pandas.DataFrame` with columns (minimally) ['mol','atag','btag',feature]
    :param feature: feature by which to sort the data by
    :type feature: str (must be a column name of data)

    :return: dataframe sorted by feature
    :rtype: :class:`pandas.DataFrame`
    """

    if not isinstance(data, pd.DataFrame):
        raise TypeError('data provided needs to be a pandas DataFrame object')

    data.sort_values(by=feature, axis=0, ascending=False, inplace=True)  # sort by specified feature

    return data


def get_selection(r_group, no_picks, data):
    """Selects a specified number of either R1 or R2 groups from a provided dataframe of groups.
    e.g. r_group_decomp.csv

    :param r_group: either 'R1' or 'R2' etc.
    :type r_group: string
    :param no_picks: number of R groups to retrieve
    :type no_picks: int
    :param data: Pandas dataframe containing all the R1 and R2 groups
    :type data: :class:`pandas.Dataframe`

    :return: [[R group identifiers], [R group SMILES]]
    :rtype: nested list
    """
    RDLogger.DisableLog('rdApp.*')  # don't print warning message for hydrogen removal

    r_groups = data[r_group.upper()]  # subset the data to only the given r group (either 'R1' or 'R2')
    mols = [Chem.MolFromSmiles(mol) for mol in r_groups]
    fps = [Chem.AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in mols]  # make fingerprints
    picker = rdSimDivPickers.MaxMinPicker() # use maxmin alg to make diverse selection
    picks = list(picker.LazyBitVectorPick(fps, len(fps), no_picks))
    r_ind = data.columns.get_loc(r_group)
    return list(data.iloc[:, r_ind-1][picks]), list(data.iloc[:, r_ind][picks])
