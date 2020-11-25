from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem import AllChem
from rdkit import RDLogger
import pandas as pd

def read_mols(filename):
    """ Reads the original csv and extracts the molecules, tags and pIC50 values

    :param filename: path to the .csv file containing the data e.g. 2010 ACS MedChemLett SI.csv
    :type filename: str

    :return: :class:`pandas.DataFrame` with columns ['mol','atag','btag','pic50']
    :rtype: :class:`pandas.DataFrame`
    """
    data = pd.read_csv(filename, delimiter=',')
    data.columns = map(str.lower, data.columns)

    # Get the relevant columns from the data
    atags = data['atag']
    btags = data['btag']
    mols = [Chem.MolFromSmiles(x) for x in data['smiles']]  # get all the mols from the smiles
    pic50 = data[data.filter(like='pic50').columns[0]]  # get the column with the pIC50 information
    mols_df = pd.DataFrame({'mol':mols, 'atag':atags, 'btag':btags, 'pic50':pic50})  # create df with extracted info

    return mols_df

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

# x = read_mols('../../data/2010 ACS MedChemLett SI.csv')
# print(x)
# x = rank_mols(x, 'pic50')
# print(x)

def get_r(data, scaffold):
    """Copy the get_r function from Get_r_groups.py, it should decompose the R groups
    from the molecules

    :param data: a dataframe containing, minimally, a list of smiles and atags and btags
    """
    scaffold = Chem.MolFromSmiles(scaffold)
    mols = data['mol']
    groups, _ = rdRGD.RGroupDecompose([scaffold], mols, asSmiles=True)
    groups_frame = pd.DataFrame(groups)
    R_cols = [col for col in groups_frame.columns if 'R' in col]
    for col in R_cols:
        num = int(col[1])
        letter = chr(ord('a') + num - 1)
        tag = letter + 'tag'
        groups_frame.insert(
            loc=groups_frame.columns.get_loc(col),
            column=tag,
            value=groups_frame[col].factorize()[0] + 1
        )


def get_selection():
    """Copy the get_selection file to return R groups that can go into BasicGameLoop.py script
    """


"""
Would be good to keep the output of the get_selection function the same so
that it can be used by Olivia's script.
"""