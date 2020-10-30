import pandas as pd
import rdkit
from combine import MolChoose
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw

R1list = ["A01", "A11", "A21"]
R2list = ["B14", "B11", "B21"]


class game:
    def __init__(self, rgroupslist=["R1", "R2"], decompFile="/home/sabsr3/Softwear Engineering Module/MainProject/localrepo/drug_discovery_game/data/r_groups_pIC50.csv"):
        print("New Game Started")
        print("to exit the game at enter \"Exit\" ")
        self.decompFile = decompFile
        self.Options = []
        # self.Options.append(self.getOptions("R1"))
        # self.Options.append(self.getOptions(R2list, "R2"))
        self.rgroupslist = rgroupslist
        self.currrentchoice = ["", ""]
        self.currentoptions = ["", ""]
        self.play()

    def play(self,):
        self.Exit = False
        while self.Exit is False:
            for group in range(len(self.rgroupslist)):
                self.currentoptions[group] = self.getOptions(self.rgroupslist[group])
            self.makeChoice(self.currentoptions[group])
            self.giveFeedback()
    

    def getOptions(self, GroupSelect="R1"):
        """

        """

        RList = self.steffunciton(GroupSelect)
        try:
            decomp = pd.read_csv(self.decompFile)
        except FileNotFoundError:
            print("file specified ("+str(self.decompFile)+") appers to not exist. the file must exist")
            self.decompFile = input()
        decomp = pd.read_csv(self.decompFile)
        rOptions = self.getSets(decomp, RList, GroupSelect)
        img = Draw.MolsToGridImage(rOptions["Mol"][:], molsPerRow=len(rOptions["Mol"]), subImgSize=(400, 400), legends=list(rOptions["tag"]))
        img.save('/home/sabsr3/Softwear Engineering Module/MainProject/localrepo/drug_discovery_game/Images/test'+str(GroupSelect)+'image.png')
        return rOptions

    def getSets(self, decomp, GroupList, GroupSelect="R1"):
        """function that when given a list of rgoup names will return a pandas dataframe with the corrosponding rgorup data from the file containing data on the rgroups.


        """
        rOptions = pd.DataFrame(columns=["tag", "Smiles"])  #setup dataframe
        if GroupSelect == "R1": #if reading first R group
            colName = "Atag"
        elif GroupSelect == "R2": # if reading second R group
            colName = "Btag"

        for option in range(len(GroupList)):
            commonOtherR = decomp.loc[decomp[colName] == GroupList[option]]
            commonOtherR = commonOtherR.iloc[0]
            # r1mol = Chem.MolFromSmiles(commonOtherR["R1"])
            # r1mol = Chem.rdmolops.DeleteSubstructs(r1mol, Chem.MolFromSmiles("[H][*]"))
            # commonOtherR["R1"] = Chem.MolToSmiles(r1mol)
            rAtributes = [commonOtherR[colName], commonOtherR[GroupSelect]]
            rOptions.loc[option] = rAtributes
        rdkit.Chem.PandasTools.AddMoleculeColumnToFrame(rOptions, smilesCol='Smiles', molCol='Mol')
        return rOptions

    def makeChoice(self,options):
        print("options")
        print(options)
        validChoice = False
        while validChoice == False:
            for i in range(len(self.rgroupslist)):
                Rchoice = input("Choice for "+str(self.rgroupslist[i])+" (number only): ")
                if ("Exit" in Rchoice) or ("exit" in Rchoice) or (Rchoice == "E") or (Rchoice == "e"):
                    self.Exit = True
                    return
                else:
                    self.currrentchoice[i] = Rchoice

            print("you have selected: " + str(self.currrentchoice))
            self.chosenmolecule = "a" # MolChoose(self.currrentchoice[0], self.currrentchoice[1])
            if str(type(self.chosenmolecule)) == "<class 'pandas.core.frame.DataFrame'>":
                validChoice = True
            else:
                print("That choice was invalid:")
                print(self.chosenmolecule)
                retry = input("do you want to reselect again? [\033[4mY\033[0mes/no]")
                if ("Exit" in retry) or ("exit" in retry) or ("No" in retry) or ("no" in retry) or (retry == "E") or (retry == "e"):
                    self.Exit = True
                    return
        return self.chosenmolecule
        


    def giveFeedback(self):
        # placeholder
        print("giveFeedback() is a place holder")

    def steffunciton(self, group):
        if group == "R1":
            return ["A01", "A11", "A21"]
        elif group == "R2":
            return ["B14", "B11", "B21"]


game()
