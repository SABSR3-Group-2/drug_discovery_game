"""Tests that CSV data is imported correctly"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import rdRGroupDecomposition as rdRGD
import pandas as pd
from io import StringIO

def test_csv_read
