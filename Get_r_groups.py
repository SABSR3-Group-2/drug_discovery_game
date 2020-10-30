from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition as rdRGD
import pandas as pd


class Get_r_groups:
    """
    Returns a dataframe containing R groups from a .csv file of molecules.

    Parameters
    :param filename: The name of the csv file
    :type filename: str
    :param smiles: The smiles representation of the molecular scaffold
    :type smiles: str
    """
    def __init__(self, filename, smiles):
        self.filename = filename
        self.smiles = smiles

    def get_r(self):
        try:
            data = pd.read_csv(self.filename, delimiter=',')
            data.columns = map(str.lower, data.columns)
            scaffold = Chem.MolFromSmiles(self.smiles)
            mols = [Chem.MolFromSmiles(x) for x in data['smiles']]
            for col in data.columns:
                if 'pic50' in col:
                    data.rename(columns={col: 'pic50'}, inplace=True)
            groups, _ = rdRGD.RGroupDecompose([scaffold], mols, asSmiles=True)
            groups_frame = pd.DataFrame(groups)
            for new_col, index in zip(['atag', 'btag', 'pic50'], [1, 3, 5]):
                groups_frame.insert(index, new_col, data[new_col])
        except FileNotFoundError:
            print("File specified " + self.filename + " does not exist.")
        except KeyError:
            raise RuntimeError(self.filename + " mising correct column names.")
        return groups_frame
