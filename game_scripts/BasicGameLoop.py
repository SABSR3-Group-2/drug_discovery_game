import pandas as pd
import rdkit
from combine import MolChoose
from Get_r_groups import Get_r_groups
from r_groups_selection import get_selection
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from r_groups_selection import get_selection
#import pillow
import os
_ = os.system('clear') 

R1list = ["A01", "A11", "A21"]
R2list = ["B14", "B11", "B21"]


class game:
    def __init__(self, rgroupslist=["R1", "R2"], decompFile="/home/sabsr3/Softwear Engineering Module/MainProject/drug_discovery_game/data/r_group_decomp.csv"):
        self.democheat()
        print("--------------------\nNew Game Started\n--------------------")
        print("to exit the game at enter \"Exit\"")
        self.decompFile = decompFile
        self.Options = []
        # self.Options.append(self.getOptions("R1"))
        # self.Options.append(self.getOptions(R2list, "R2"))
        self.rgroupslist = rgroupslist
        self.currrentchoice = ["",""]
        self.currentoptions = ["", ""]
        self.play()

    def play(self,termlimit=5):
        self.__scores = []
        self.choiceHistory = []
        self.__loopnumber = 1
        self.Exit = False
        while self.Exit is False:
            print("\n\nRound "+str(self.__loopnumber)+" of "+str(termlimit))
            self.currrentchoice = list([])
            choice = ["",""]
            for group in range(len(self.rgroupslist)):
                self.currentoptions[group] = self.getOptions(self.rgroupslist[group])

            self.goodcombo = False
            while self.goodcombo == False:
                for group in range(len(self.rgroupslist)):
                    choice[group] = (self.makeChoice(self.currentoptions[group],self.rgroupslist[group]))
                    if self.Exit == True:
                        break
                if self.Exit == True:
                    break
                self.currrentchoice = choice
                self.molCheck()
                self.chosenmolecule = MolChoose(self.currrentchoice[0], self.currrentchoice[1])
            if self.Exit == True:
                break
            self.choiceHistory.append(self.currrentchoice)
            self.giveFeedback()
            self.__loopnumber += 1
            if self.__loopnumber > termlimit:
                self.Exit = True
                self.endgame()

            
    

    def getOptions(self, GroupSelect="R1"):
        """

        """
        #RList = self.get_selection(GroupSelect) #add self. if local
        try:
            decomp = pd.read_csv(self.decompFile)
        except FileNotFoundError:
            print("file specified ("+str(self.decompFile)+") appers to not exist. the file must exist")
            self.decompFile = input()
        decomp = pd.read_csv(self.decompFile)

        self.Rlist = [get_selection("R1", 5, decomp)[0], get_selection("R2", 5, decomp)[0]] 
        rOptions = self.getSets(decomp, self.Rlist, GroupSelect)
        img = Draw.MolsToGridImage(rOptions["Mol"][:], molsPerRow=len(rOptions["Mol"]), legends=list(rOptions["tag"]))# , subImgSize=(400, 400))
        #img.show()
        img.save('Images/test'+str(GroupSelect)+'image.png')
        #import os
        #os.system('Images/test'+str(GroupSelect)+'image.png')
        return rOptions

    def getSets(self, decomp, GroupList, GroupSelect="R1"):
        """function that when given a list of rgoup names will return a pandas dataframe with the corrosponding rgorup data from the file containing data on the rgroups.


        """
        rOptions = pd.DataFrame(columns=["tag", "Smiles"])  #setup dataframe
        if GroupSelect == "R1": #if reading first R group
            colName = "atag"
            groupnum = 0
        elif GroupSelect == "R2": # if reading second R group
            colName = "btag"
            groupnum = 1

        for option in range(len(GroupList[groupnum])):
            # print(GroupList)
            # print("colls names")
            # print(decomp[colName])
            # print("group list options")
            # print(GroupList[groupnum][option])
            # print("full thing")
            # print(decomp.loc[decomp[colName] == GroupList[groupnum][option]])
            # print("\n")
            commonOtherR = decomp.loc[decomp[colName] == GroupList[groupnum][option]]
            commonOtherR = commonOtherR.iloc[0]
            # r1mol = Chem.MolFromSmiles(commonOtherR["R1"])
            # r1mol = Chem.rdmolops.DeleteSubstructs(r1mol, Chem.MolFromSmiles("[H][*]"))
            # commonOtherR["R1"] = Chem.MolToSmiles(r1mol)
            rAtributes = [commonOtherR[colName], commonOtherR[GroupSelect]]
            rOptions.loc[option] = rAtributes
        rdkit.Chem.PandasTools.AddMoleculeColumnToFrame(rOptions, smilesCol='Smiles', molCol='Mol')
        return rOptions

    def makeChoice(self,options,groupname):
        oppdata = pd.DataFrame(options)
        validChoice = False
        while validChoice == False:
            Rchoice = input("> Choice for "+str(groupname)+" " + str(list(oppdata.loc[:,'tag'])) +": ")                          #add clearer options
            if ("Exit" in Rchoice) or ("exit" in Rchoice) or (Rchoice == "E") or (Rchoice == "e"):
                self.Exit = True
                return
            else:
                self.currrentchoice = Rchoice

            print("- you have selected: " + str(self.currrentchoice))
            self.chosenmolecule = MolChoose(self.currrentchoice[0], self.currrentchoice[1])
            if str(type(self.chosenmolecule)) == "<class 'pandas.core.frame.DataFrame'>":
                validChoice = True
            else:
                print("That choice was invalid:")
                print(self.chosenmolecule)
                retry = input("do you want to reselect again? [\033[4mY\033[0mes/no]:")
                if ("Exit" in retry) or ("exit" in retry) or ("No" in retry) or ("no" in retry) or (retry == "E") or (retry == "e"):
                     self.Exit = True
                     return 
        return Rchoice

    def molCheck(self):
        checkcombo = True
        for choice in self.currrentchoice:
            goodchoice =False
            for set in self.currentoptions:
                if choice in list(set.loc[:,'tag']):
                    goodchoice = True
            if goodchoice == False:
                print("\n\n######################## Oh, that looks wrong ########################\n\""+str(choice)+"\" was not found in the option set, please reselect your R Groups\n######################################################################\n")
                checkcombo = False
        self.goodcombo = checkcombo
        return checkcombo
        

    def UpdateCurrentImage(self, Smiles):
        print("- Updatinging current")
        print(str(Smiles[1]))
        currentmol = rdkit.Chem.MolFromSmiles(str(Smiles))
        print(currentmol)
        img = Draw.MolsToGridImage(currentmol,  molsPerRow=1, subImgSize=(400, 400))
        img.save('Images/CurrentImage.png')


    def giveFeedback(self):
        # placeholder
        #print("- giveFeedback() is a place holder")
        #print(self.chosenmolecule)
        pIC50 = self.chosenmolecule.loc[:,"pIC50"].to_string(index=False)
        print("\n- that molecule's pIC50 is "+str(pIC50))
        self.__scores.append(pIC50)
        #print((self.chosenmolecule))
        #currentstring = self.chosenmolecule.loc["Core"]
        #print(currentstring)
        #self.UpdateCurrentImage(currentstring) #.values[0]

    # def get_selection(self, group):
    #     #placeholder for full get_selection fucntion
    #     if group == "R1":
    #         return ["A02", "A03", "A04"]                                          #add smiles
    #     elif group == "R2":
    #         return ["B27", "B28", "B29"]


    def endgame(self):
        print("Thank you for playing")
        print(self.__scores)
        self.plotscores()

    def plotscores(self):
        import numpy as np
        import matplotlib.pyplot as plt
        print(self.choiceHistory)
        numericscorelist = []
        fig, ax = plt.subplots()
        for scorenum in range(len(self.__scores)):
            currentscore = self.__scores[scorenum]
            if currentscore == " Not Made":
                colour = "grey"
                currentscore = 0.0
                labelText = "Not Made"
            elif  currentscore == " Assay Failed":
                colour = "red"
                currentscore = 0.0
                labelText = "Assay Failed"
            elif currentscore == " Inactive":
                colour = "blue"
                currentscore = 0.0
                labelText = "Inactive"

            else:
                currentscore = float(str(currentscore).replace(" ",''))
                colour = "limegreen"

                labelText = "pIC50"
            print(str(scorenum)+" "+str(currentscore))
            ax.scatter(float(scorenum), float(currentscore), c=colour, label=labelText)
            ax.annotate(str(self.choiceHistory[scorenum]).replace("\'",''), (scorenum, currentscore))
            numericscorelist.append(float(currentscore))

        #ax.legend(False)
        ax.set_xlim(0,len(numericscorelist))
        ax.set_ylim(0,max(numericscorelist))

        plt.title('Chosen Molecule pIC50')
        plt.xlabel("Round")
        plt.ylabel("pIC50")
        plt.show()
        plt.savefig('Images/CurrentImage.png')


    def democheat(self):
        from rdkit import Chem                                                                                                                                                                   
        from rdkit import RDLogger                                                                                                                                                          
        RDLogger.DisableLog('rdApp.*') 

game()
