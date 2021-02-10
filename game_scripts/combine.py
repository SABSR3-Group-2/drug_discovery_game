import pandas as pd
import os.path


def MolChoose(R1,R2,DataSource="data/r_group_decomp.csv"):
    """a function to return the source data on the molecule specified by the provided R1 and R2 values
    


    :param R1: Tag specifying R1 group
    :type R1: String 
    :param R2: Tag specifying R2 group
    :type R2: String
    :param DataSource: file path to file providing full molecule smile strings, defaults to data/r_groups_pIC50.csv
    :type DataSource: String
    ...
    :raises FileNotFoundError: the file listing the data as named in DataSource is not able to be opened
    :raises RuntimeError: if the R1 column (as named in R1col) is not present in the dataset
    :raises RuntimeError: if the R2 column (as named in R2col) is not present in the dataset
    ...
    :return: selectedrow - returns data listed in original document relevent to the molecule combining R1 and R2
    :rtype: Pandas row

    """
    R1col = 'atag'
    R2col = 'btag'
    try:
        originalset = pd.read_csv(DataSource)
    except FileNotFoundError as e:   
        print("The specified file ("+str(DataSource)+") cannot be found")
        raise e
    try:
        selectedrows = originalset.loc[originalset[R1col] == R1]
    except KeyError as k:
        print('No \""+str(R1col)+"\" in \""+str(DataSource)+"\"')
        raise k
    try:
         selectedrow = selectedrows.loc[selectedrows[R2col] == R2]
    except KeyError as k:
        print('No \""+str(R2col)+"\" in \""+str(DataSource)+"\"')
        raise k
    return selectedrow
