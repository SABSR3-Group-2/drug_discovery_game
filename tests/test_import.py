import unittest
import drug_discovery_game
from rdkit.Chem import PandasTools

class ImportTest(unittest.TestCase):
    """
    Tests that the functions for importing chemical information from input files is working
    """
    def test_import_format(self):
