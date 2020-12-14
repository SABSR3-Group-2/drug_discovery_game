from rdkit import Chem
from rdkit.Chem import FilterCatalog
import pandas as pd
from rdkit.Chem.FilterCatalog import *

# function to filter a list of molecules
def pains_filter(mols):
    """Filters pan-assay interference compounds (PAINS) using RDKit's
    filter catalog.

    :param mols: list of RDKit molecules to filter
    :type mols: list

    :return: 
    :rtype:
    """
    pains_compounds = []
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog(params)
    for mol in mols:
        if catalog.HasMatch(mol):
            pains_compounds.append(mol)
    return pains_compounds

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