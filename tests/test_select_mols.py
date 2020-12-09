"""Tests that the select_mols.py script functions work as expected"""

import pytest
import pandas as pd
import os.path
from game_scripts import select_mols

cwd = os.path.dirname(__file__)  # get current working directory


# Test the get selection function
def test_correct_atag():
    """Tests that the correct 'tags' are returned given the specified R group. e.g. A1 would be a correct return for R1
    and B3 for R2 etc."""

    test_data = pd.read_csv(os.path.join(cwd, 'test_get_selec.csv'))
    test_selection = select_mols.get_selection('R1', 3, test_data)
    assert (tag[0] == 'A' for tag in test_selection[0])


def test_correct_number():
    """Tests that the number of R groups returned is equal to the specified number"""

    test_data = pd.read_csv(os.path.join(cwd, 'test_get_selec.csv'))
    test_selection = select_mols.get_selection('R2', 4, test_data)
    assert len(test_selection[0]) == 4


# # Test the read_mols function
# def test_read_mols():
#     """Tests that a dataframe is returned with the correct columns"""
#
#     select_mols.read_mols()
