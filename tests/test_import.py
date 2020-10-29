import unittest
from Get_r_groups import Get_r_groups
from rdkit.Chem import PandasTools
from rdkit import Chem
import pandas as pd


class ImportTest(unittest.TestCase):
    """
    Tests that the functions for importing chemical information from input files is working
    """

    data = {'Index': [1, 2, 3],
            'pIC50_MMP12': [6.5, 6.2, 5.0],
            'Smiles': ['OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3ccccc3',
                       'OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3ccc(Br)cc3',
                       'OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3cccc(c3)N(=O)=O'],
            'Final QC Purity': ['>98', '>98', '>97'],
            'Atag': ['A01', 'A01', 'A01'],
            'Btag': ['B01', 'B02', 'B03']}

    def test_ROMol(self, data=data):
        """
        Test for checking that MolFromSmiles works correctly
        :param data: Dataframe containing some information on 3 molecules
        :type data: pandas dataframe object
        :return:
        """
        df = pd.DataFrame.from_dict(data)
        mols = [Chem.MolFromSmiles(x) for x in df['Smiles']]
        test_mol = Chem.MolFromSmiles('OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3ccccc3')
        types = [type(x) for x in mols]
        test_types = [type(test_mol) for i in range(3)]
        self.assertEqual(types, test_types)

    def test_get_r(self, data=data):
        """
        Test the full get_r() function from the Get_r_groups script

        :param data:
        :return:
        """

        test_r =