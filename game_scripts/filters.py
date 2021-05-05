from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import *


def compound_filter(smiles):
    """Filters undesirable compounds using RDKit's filter catalog.
    Includes PAINS, ZINC, BRENK and NIH filters.

    :param smiles: smiles of molecule for filtering
    :type mols: list

    :return: true or false indicating if undesirable compound
    :rtype: Boolean value
    """

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    catalog = FilterCatalog(params)
    mol = Chem.MolFromSmiles(smiles)
    if catalog.HasMatch(mol):
        return False
    else:
        return True


def apply_filter(df, r_group):
    """Filters a dataframe for undesirable compounds using the filter function.

    :param df: dataframe containing molecules for filtering
    :type df: Pandas DataFrame
    :param r_group: the r_group to run the filtering on
    :type r_group: string

    :return: filtered dataframe containing passing compounds
    :rtype: PandasDataFrame
    """
    filtered_df = df[df.apply(lambda x: compound_filter(x[r_group]), axis=1)]
    return filtered_df


def compound_check(mol):
    """Checks whether a molecule is an undesirable compound and
    identifies the reason and substructure responsible.

    :param mols: molecule to check
    :type mols: RDKit molecule

    :return: warning if molecule is an undesirable compound and why
    :rtype: string
    """

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    catalog = FilterCatalog(params)
    entries = catalog.GetMatches(mol)
    results = []
    if entries:
        no = 1
        for entry in entries:
            result = f'Match {no}\nWarning: molecule failed filter\nReason: {entry.GetDescription()}'
            results.append(result)
            no += 1
        return results
    else:
        result = 'Molecule passes the filter.'
        results.append(result)
        return results

def run_filters(mol):
    """Checks whether a molecule is an undesirable compound against each filter and
    identifies the reason and substructure responsible.

    :param mols: molecule to check
    :type mols: RDKit molecule

    :return: dictionary containing filter name as key and warnings/passing as values
    :rtype: dict
    """

    pains_params = FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)

    zinc_params = FilterCatalogParams()
    zinc_params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)

    brenk_params = FilterCatalogParams()
    brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)

    nih_params = FilterCatalogParams()
    nih_params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)

    params = [pains_params, zinc_params, brenk_params, nih_params]
    filt_names = ['PAINS', 'ZINC', 'BRENK', 'NIH']

    all_warnings = {}

    for param, name in zip(params, filt_names):
        catalog = FilterCatalog(param)
        entries = catalog.GetMatches(mol)
        warnings = []
        for e in entries:
            warnings.append(e.GetDescription())
        if len(warnings) != 0:
            all_warnings[name] = warnings
        else:
            all_warnings[name] = 'passing'
    
    return all_warnings
        