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
                # rename pic50 column to just pic50
            test_frame = {'Core': mols, 'R1': mols, 'R2': mols}
            groups_frame = pd.DataFrame(test_frame, columns=['Core', 'R1', 'R2'])
            for new_col, index in zip(['atag', 'btag', 'pic50'], [1, 3, 5]):
                groups_frame.insert(index, new_col, data[new_col])
        except FileNotFoundError:
            print("File specified " + self.filename + " does not exist.")
        except KeyError:
            raise RuntimeError("No smiles column in " + self.filename)
        return groups_frame

Get_r_groups('spreadsheet2.csv', 'O=S(C1=CC=CC=C1)(NCC(O)=O)=O').get_r()
