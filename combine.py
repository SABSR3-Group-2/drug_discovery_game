import pandas as pd
import os.path


def MolChoose(R1,R2,DataSource="data/2011.Pickett-et-al.AutomatedLeadOptimizationOfMMP-12InhibitorsUsingAGeneticAlgorithm.csv"):
    """a function to return the source data on the molecule specified by the provided R1 and R2 values
    


    :param R1: Tag specifying R1 group
    :type R1: String 
    :param R2: Tag specifying R2 group
    :type R2: String
    :param DataSource: file path to file providing full molecule smile strings, defaults to data/2011.Pickett-et-al.AutomatedLeadOptimizationOfMMP-12InhibitorsUsingAGeneticAlgorithm.csv
    :type DataSource: String
    ...
    :raises FileNotFoundError: the file listing the data as named in DataSource is not able to be opened
    :raises RuntimeError: if the R1 column (as named in R1col) is not present in the dataset
    :raises RuntimeError: if the R2 column (as named in R2col) is not present in the dataset
    ...
    :return: selectedrow - returns data listed in original document relevent to the molecule combining R1 and R2
    :rtype: Pandas row

    """
    R1col = 'Atag'
    R2col = 'Btag'
    try:
        origninalset = pd.read_csv(DataSource)
    except FileNotFoundError:
        print("file specified ("+str(DataSource)+") appers to not exist. the file must exist")
    try:
        selectedrows = origninalset.loc[origninalset[R1col] == R1]
    except KeyError:
        raise RuntimeError("no \""+str(R1col)+"\" collumn in \""+str(DataSource)+"\"")
    try:
         selectedrow = selectedrows.loc[selectedrows[R2col] == R2]
    except KeyError:
        raise RuntimeError("no \""+str(R2col)+"\" collumn in \""+str(DataSource)+"\"")
    
    origninalset = pd.read_csv(DataSource)
    selectedrows = origninalset.loc[origninalset['Abtag'] == R1]
    selectedrow = selectedrows.loc[selectedrows['Btag'] == R2]
    return selectedrow

print(MolChoose("A06","B01"))