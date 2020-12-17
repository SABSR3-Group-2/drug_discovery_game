import numpy.testing as npt
import pandas as pd
from game_scripts import filters
from rdkit import Chem


def test_compound_filter():
    """Tests a molecule containing a thiocarbonyl fails the filter"""
    from game_scripts.filters import compound_filter
    pains_mol = 'O=C(NC(C1=C(C)C=C2N1C=CC=C2)=S)C3=CC=CC=C3'
    npt.assert_equal(False, compound_filter(pains_mol))


def test_apply_filter():
    """Tests that a dataframe is filtered correctly"""
    from game_scripts.filters import apply_filter
    test_df = pd.DataFrame({'mols': [
        'C=CCN1C(CCC(O)=O)=CC=C1C2=CC=C(F)C=C2',
        'OC1=C(C(O)=O)C=CC=C1O',
        'CC'
        ]})
    correct_df = pd.DataFrame({'mols': ['CC']})
    correct_df.index = correct_df.index + 2
    pd.testing.assert_frame_equal(correct_df, apply_filter(test_df, 'mols'))


def test_compound_check():
    """Tests reason for failing the filters is correctly identified"""
    from game_scripts.filters import compound_check
    pains_mol = Chem.MolFromSmiles('O=C(C1=C2C(C)OC(CC(O)=O)C1)C3=CC=CC(O)=C3C2=O')
    warning =  """
    Match  1
    Warning: molecule failed filter: reason quinone_A(370)
    PAINS filters (family A)
    """
    npt.assert_equal(print(warning), compound_check(pains_mol))


def test_compound_check_2():
    """Tests a molecule correctly passes the filters"""
    from game_scripts.filters import compound_check
    message = 'Molecule passes the filter.'
    passing_mol = Chem.MolFromSmiles('CC')
    npt.assert_equal(print(message), compound_check(passing_mol))
