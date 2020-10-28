from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition as rdRGD
import pandas as pd


class Get_r_groups:
    """
    Returns a dataframe containing R groups of a .csv file of molecules.

    Parameters
    :filename: the name of the csv file
    :filename type: string
    :smiles: the smiles representation of the molecular scaffold
    :smiles type: string
    """
    def __init__(self, filename, smiles):
        self.filename = filename
        self.smiles = smiles

    def get_r(self):
        # will tidy up exception handling
        try:
            data = pd.read_csv(self.filename, delimiter=',')
        except FileNotFoundError:
            print("File specified " + self.filename + "does not exist.")
        scaffold = Chem.MolFromSmiles(self.smiles)
        try:
            mols = [Chem.MolFromSmiles(x) for x in data['Smiles']]
        except KeyError:
            raise RuntimeError("No 'Smiles' column in " + self.filename)
        groups, _ = rdRGD.RGroupDecompose([scaffold], mols, asSmiles=True)
        groups_frame = pd.DataFrame(groups)
        try:
            groups_frame.insert(1, 'Atag', data['Atag'])
            groups_frame.insert(1, 'Btag', data['Btag'])
            groups_frame.insert(1, 'pIC50', data['pIC50_MMP12'])
        except KeyError:
            raise RuntimeError("R group and pIC50 columns missing from file")
        return groups_frame
