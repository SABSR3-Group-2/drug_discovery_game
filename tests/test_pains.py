import numpy.testing as npt
import pandas as pd
from game_scripts import pains
from rdkit import Chem


def test_pains_filter():
    """Tests a molecule containing a PAINS substructure fails the filter"""
    from game_scripts.pains import pains_filter
    pains_mol = 'O=C(NC(C1=C(C)C=C2N1C=CC=C2)=S)C3=CC=CC=C3'
    npt.assert_equal(False, pains_filter(pains_mol))


def test_apply_filter():
    """Tests that a dataframe is filtered correctly"""
    from game_scripts.pains import apply_filter
    test_df = pd.DataFrame({'mols': [
        'C=CCN1C(CCC(O)=O)=CC=C1C2=CC=C(F)C=C2',
        'O=C(C1=C2C(C)OC(CC(O)=O)C1)C3=CC=CC(O)=C3C2=O',
        'CC'
        ]})
    correct_df = pd.DataFrame({'mols': ['CC']})
    correct_df.index = correct_df.index + 2
    pd.testing.assert_frame_equal(correct_df, apply_filter(test_df, 'mols'))


def test_pains_check():
    """Tests reason for failing the PAINS filter is correctly identified"""
    from game_scripts.pains import pains_check
    warning = 'Warning: molecule failed filter: reason catechol_A(92)'
    pains_mol = Chem.MolFromSmiles('OC1=C(C(O)=O)C=CC=C1O')
    npt.assert_equal(print(warning), pains_check(pains_mol))


def test_pains_check_2():
    """Tests a molecule correctly passes PAINS filter"""
    from game_scripts.pains import pains_check
    message = 'Molecule passes the PAINS filter.'
    passing_mol = Chem.MolFromSmiles('CC')
    npt.assert_equal(print(message), pains_check(passing_mol))
