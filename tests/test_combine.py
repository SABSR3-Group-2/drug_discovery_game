"""Tests that the combine.py script functions as expected"""

import pytest
import unittest
import pandas as pd
import os.path
import combine

def test_molchoose_nofile():
    """Tests that the FileNotFoundError works as expected"""
    from combine import MolChoose 

    with pytest.raises(FileNotFoundError):
        MolChoose('C', 'c', DataSource ='data/no_such_file.csv')

test_molchoose_nofile()
