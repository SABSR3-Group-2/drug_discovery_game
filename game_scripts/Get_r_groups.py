from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem import AllChem
from rdkit import RDLogger
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

    def remove_h(self, r1):
        RDLogger.DisableLog('rdApp.*')
        if r1 == '[H][*:1].[H][*:1]':
            return '[H][*:1]'
        else:
            mol = Chem.MolFromSmiles(r1)
            substruct = Chem.MolFromSmiles('[H][*:1]')
            mol = AllChem.DeleteSubstructs(mol, substruct)
            Chem.SanitizeMol(mol)
            return Chem.MolToSmiles(mol)

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
                groups_frame[tag] = groups_frame[tag].apply("{:02d}".format)
                groups_frame[tag] = letter.upper() + groups_frame[tag].astype(str)
                groups_frame[col] = groups_frame[col].apply(self.remove_h)
            groups_frame['pic50'] = data['pic50']
        except FileNotFoundError:
            print("File specified " + self.filename + " does not exist.")
        except KeyError:
            raise RuntimeError(self.filename + " missing correct col names.")
        return groups_frame
