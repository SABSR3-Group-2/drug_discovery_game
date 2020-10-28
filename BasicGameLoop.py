import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw

R1list = ["A01","A11","A21","A31","A41"]
R2list = ["B01","B11","B21","B31","B41"]

brokendownfile = "localrepo/drug_discovery_game/data/r_groups_pIC50.csv"
R1col = "Atag"
try:
    brokendown = pd.read_csv(brokendownfile)
except FileNotFoundError:
    print("file specified ("+str(brokendownfile)+") appers to not exist. the file must exist")
brokendown = pd.read_csv(brokendownfile)

R1options = pd.DataFrame(columns=["tag","Smiles"])
i=0
for Group in R1list:
    commonR1 = brokendown.loc[brokendown[R1col] == Group]
    commonR1 = commonR1.iloc[0]
    R1atributes = [commonR1[R1col],commonR1["R1"]] #pd.DataFrame((commonR1[R1col],commonR1["R1"]))
    R1options.loc[i] = R1atributes
    i+=1
print(R1options)
PandasTools.AddMoleculeColumnToFrame(R1options, smilesCol='Smiles', molCol='Mol')

print(R1options)
img=Draw.MolsToGridImage(R1options["Mol"][:],molsPerRow=1,subImgSize=(200,200),legends=["1","2","3","4","5"])#list(R1options["tag"])) 
img.save('Images/testR1image.png') 