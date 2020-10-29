import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw

R1list = ["A01","A11","A21","A31","A41"]
R2list = ["B01","B11","B21","B31","B41"]

decompFile = "/home/sabsr3/Softwear Engineering Module/MainProject/localrepo/drug_discovery_game/data/r_groups_pIC50.csv"
R1col = "Atag"
try:
    decomp = pd.read_csv(decompFile)
except FileNotFoundError:
    print("file specified ("+str(decompFile)+") appers to not exist. the file must exist")
decomp = pd.read_csv(decompFile)

R1options = pd.DataFrame(columns=["tag","Smiles"])
i=0
for Group in R1list:
    commonR1 = decomp.loc[decomp[R1col] == Group]
    commonR1 = commonR1.iloc[0]
    R1atributes = [commonR1[R1col],commonR1["R1"]] #pd.DataFrame((commonR1[R1col],commonR1["R1"]))
    R1options.loc[i] = R1atributes
    i+=1
PandasTools.AddMoleculeColumnToFrame(R1options, smilesCol='Smiles', molCol='Mol')

img=Draw.MolsToGridImage(R1options["Mol"][:],molsPerRow=1,subImgSize=(200,200),legends=["1","2","3","4","5"])#list(R1options["tag"])) 
img.save('/home/sabsr3/Softwear Engineering Module/MainProject/localrepo/drug_discovery_game/Images/testR1image.png') 


def getOptions(decompFile = "/home/sabsr3/Softwear Engineering Module/MainProject/localrepo/drug_discovery_game/data/r_groups_pIC50.csv"):
    R1list = ["A01","A11","A21","A31","A41"]
    R2list = ["B01","B11","B21","B31","B41"]
    try:
        decomp = pd.read_csv(decompFile)
    except FileNotFoundError:
        print("file specified ("+str(decompFile)+") appers to not exist. the file must exist")
        decompFile = input()
    decomp = pd.read_csv(decompFile)
    getSets(decomp,R2list,"R2")


def getSets(decomp, GroupList, GroupSelect="R1"):
    R1options = pd.DataFrame(columns=["tag","Smiles"]) #setup dataframe
    if GroupSelect == "R1": #if reading first R group
        colName = "Atag"
    elif GroupSelect == "R2": # if reading second R group
        colName = "Btag"

    for option in range(len(GroupList)):
        commonOtherR = decomp.loc[decomp[colName] == GroupList[option]]
        print("commonOtherR")
        print(commonOtherR)
        commonOtherR = commonOtherR.iloc[0]
        R1atributes = [commonOtherR[colName],commonOtherR[GroupSelect]]
        R1options.loc[option] = R1atributes
    PandasTools.AddMoleculeColumnToFrame(R1options, smilesCol='Smiles', molCol='Mol')


getOptions()