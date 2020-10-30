import unittest
from Get_r_groups import Get_r_groups
from rdkit.Chem import PandasTools
from rdkit import Chem
import pandas as pd
from pandas.testing import assert_frame_equal


class ImportTest(unittest.TestCase):
    """
    Tests that the functions for importing chemical information from input files is working
    """

    def test_ROMol(self):
        """
        Test for checking that MolFromSmiles works correctly
        :param data: Dataframe containing some information on 3 molecules
        :type data: pandas dataframe object
        :return:
        """
        data = {'Index': [1, 2, 3],
                'pIC50_MMP12': [6.5, 6.2, 5.0],
                'Smiles': ['OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3ccccc3',
                           'OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3ccc(Br)cc3',
                           'OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3cccc(c3)N(=O)=O'],
                'Final QC Purity': ['>98', '>98', '>97'],
                'Atag': ['A01', 'A01', 'A01'],
                'Btag': ['B01', 'B02', 'B03']}

        df = pd.DataFrame.from_dict(data)
        mols = [Chem.MolFromSmiles(x) for x in df['Smiles']]
        test_mol = Chem.MolFromSmiles('OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3ccccc3')
        types = [type(x) for x in mols]
        test_types = [type(test_mol) for i in range(3)]
        self.assertEqual(types, test_types)

    def test_get_r(self):
        """
        Tests that the :class:`Get_r_groups` output is correctly being made by the `get_r()` function.

        """

        output_csv = Get_r_groups('test_input.csv', 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]').get_r()
        test_output = pd.DataFrame.from_dict({'Core': ['O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]',
                                                       'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]',
                                                       'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]',
                                                       'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]',
                                                       'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]'],
                                              'pIC50': ['6.5', '6.8', 'Assay Failed', '7.1', '5.1'],
                                              'Btag': ['B01', 'B02', 'B03', 'B04', 'B05'],
                                              'Atag': ['A01', 'A01', 'A01', 'A01', 'A01'],
                                              'R1': ['Oc1ccc(C[*:1])cc1',
                                                     'Oc1ccc(C[*:1])cc1',
                                                     'Oc1ccc(C[*:1])cc1',
                                                     'Oc1ccc(C[*:1])cc1',
                                                     'Oc1ccc(C[*:1])cc1'],
                                              'R2': ['c1ccc([*:2])cc1',
                                                     'Brc1ccc([*:2])cc1',
                                                     'O=[N+]([O-])c1cccc([*:2])c1',
                                                     'Cc1ccc([*:2])cc1',
                                                     'c1ccc2c(c1)oc1c([*:2])cccc12']
                                              })
        assert_frame_equal(output_csv, test_output)
