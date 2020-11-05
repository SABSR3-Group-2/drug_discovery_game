"""Tests that the combine.py script functions as expected"""

import pytest
import pandas as pd
import os.path
from game_scripts import combine
from io import StringIO

def test_molchoose_nofile():
    """Tests that the FileNotFoundError works as expected"""

    with pytest.raises(FileNotFoundError):  #tests that an error is raised if a non-existent file is passed to it
        combine.MolChoose('C', 'c', DataSource ='data/no_such_file.csv')

def test_molchoose_correct():
    """Tests that the correct data is returned from a trial csv document"""

    test = [[1,'A01B01','A01','B01','6.5','OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3ccccc3',7,12,'>98','Gen-5']]
    test_frame = pd.DataFrame(test, columns=['Index','Tag','Atag','Btag','pIC50_MMP12','Smiles','A_SortMax','B_SortMax','Final QC Purity','Generation-No'])

    pd.testing.assert_frame_equal(combine.MolChoose('A01', 'B01', DataSource='tests/test_input.csv'), test_frame)
 
test_molchoose_nofile()
test_molchoose_correct()
