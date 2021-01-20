"""Tests for getting descriptors for molecules and ranking them"""

import pytest
import os.path
from game_scripts import descriptors

cwd = os.path.dirname(__file__)  # get current working directory

def test_descriptors():
    """Can get all specified descriptors for single moiety and return as dict"""

    mol = 'Oc1ccc(C[*:1])cc1'
    true_case = {'mol': 'Oc1ccc(C[*:1])cc1',
                 'MW': '107.0497',
                 'logP': 1.4415,
                 'TPSA': 20.23,
                 'HA': 8,
                 'h_acc': 1,
                 'h_don': 1,
                 'rings': 1
                 }
    test_case = descriptors.get_descriptors(mol=mol)
    assert true_case == test_case