from rdkit import Chem
from rdkit.Chem import FilterCatalog
import pandas as pd
from rdkit.Chem.FilterCatalog import *

# function to filter a list of molecules
def pains_filter(smiles):
    """Filters pan-assay interference compounds (PAINS) using RDKit's
    filter catalog.

    :param smiles: 
    :type mols: list

    :return: true or false indicating if a PAINS compound
    :rtype: Boolean value
    """

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog(params)
    mol = Chem.MolFromSmiles(smiles)
    if catalog.HasMatch(mol):
        return False
    else:
        return True

def apply_filter(df, r_group):
    """Filters a dataframe for PAINS compounds using the pains_filter function.

    :param df: dataframe containing molecules for filtering
    :type df: Pandas DataFrame
    :param r_group: the r_group to run the filtering on
    :type r_group: string

    :return: filtered dataframe containing non-PAINS compounds
    :rtype: PandasDataFrame
    """
    filtered_df = df[df.apply(lambda x : pains_filter(x[r_group]),axis=1)]
    return filtered_df

# function to run on a chosen molecule and tell you if/why the molecule is a pains compound
def pains_check(mol):
    """Checks whether a chosen molecule is a PAINS compound and
    identifies the substructure responsible.

    :param mols: molecule to check
    :type mols: RDKit molecule

    :return: warning if molecule is a PAINS compounds and why
    :rtype: string
    """
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog(params)
    entry = catalog.GetFirstMatch(mol)
    if entry:
        print('Warning: molecule failed filter: reason %s'%(
            entry.GetDescription()))
    else:
        print('Molecule passes the PAINS filter.')
