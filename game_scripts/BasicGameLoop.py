import pandas as pd
import rdkit
from game_scripts.combine import MolChoose
import os
from game_scripts.select_mols import get_selection
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
# import pillow
# from Get_r_groups import Get_r_groups

# clear terminal
_ = os.system('clear')

assays = {
    'pIC50': {'cost': 70, 'duration': 0.5},
    'clearance_mouse': {'cost': 7000, 'duration': 3},
    'clearance_human': {'cost': 9000, 'duration': 3.5},
    'logd': {'cost': 1000, 'duration': 1.5},
    'pampa': {'cost': 700, 'duration': 1}
    }
# check pic50 assay duration
class game:
    """[Summary]

    :param [ParamName]: [ParamDescription], defaults to [DefaultParamVal]
    :type [ParamName]: [ParamType](, optional)
    ...
    :raises [ErrorType]: [ErrorDescription]
    ...
    :return: [ReturnDescription]
    :rtype: [ReturnType]
    """
    def __init__(self, rgroupslist=["R1", "R2"], decompFile="data/r_group_decomp.csv"):
        """[Summary]

        :param [ParamName]: [ParamDescription], defaults to [DefaultParamVal]
        :type [ParamName]: [ParamType](, optional)
        ...
        :raises [ErrorType]: [ErrorDescription]
        ...
        :return: [ReturnDescription]
        :rtype: [ReturnType]
        """
        self.democheat()
        print("--------------------\nNew Game Started\n--------------------")
        print("to exit the game type \"Exit\"")
        self.decompFile = decompFile
        self.Options = []
        # self.Options.append(self.getOptions("R1"))
        # self.Options.append(self.getOptions(R2list, "R2"))
        self.rgroupslist = rgroupslist
        self.currentchoice = ["",""]
        self.currentoptions = ["", ""]
        self.assay = ""
        self.play()

    def play(self):
         """[Summary]

        :param [ParamName]: [ParamDescription], defaults to [DefaultParamVal]
        :type [ParamName]: [ParamType](, optional)
        ...
        :raises [ErrorType]: [ErrorDescription]
        ...
        :return: [ReturnDescription]
        :rtype: [ReturnType]
        """
        self.__scores = []
        self.assay_cost = 0 # updated depending on the choice of assay
        self.start_balance = 230 # arbitrary starting balance
        self.current_balance = self.start_balance
        self.choiceHistory = []
        self.__loopnumber = 1
        self.Exit = False
        self.UpdateCurrentImage()
        while self.Exit is False:
            print("\n\nRound "+str(self.__loopnumber))
            self.giveStartBalance()
            self.currentchoice = list([])
            choice = ["",""]
            for group in range(len(self.rgroupslist)):
                self.currentoptions[group] = self.getOptions(self.rgroupslist[group])

            self.goodcombo = False
            while self.goodcombo is False:
                for group in range(len(self.rgroupslist)):
                    choice[group] = (self.makeChoice(self.currentoptions[group], self.rgroupslist[group]))
                    if self.Exit is True:
                        break
                if self.Exit is True:
                    break
                self.currentchoice = choice
                self.molCheck()
                self.chosenmolecule = MolChoose(self.currentchoice[0], self.currentchoice[1])
            if self.Exit is True:
                break
            self.choiceHistory.append(self.currentchoice)
            self.chooseAssay()
            self.__loopnumber += 1

    def getOptions(self, GroupSelect="R1"):
        """
        """
        # RList = self.get_selection(GroupSelect) #add self. if local
        try:
            decomp = pd.read_csv(self.decompFile)
        except FileNotFoundError:
            print("file specified ("+str(self.decompFile)+") appears to not exist. the file must exist")
            self.decompFile = input()
        decomp = pd.read_csv(self.decompFile)

        self.Rlist = [get_selection("R1", 3, decomp)[0], get_selection("R2", 3, decomp)[0]]
        rOptions = self.getSets(decomp, self.Rlist, GroupSelect)
        img = Draw.MolsToGridImage(rOptions["Mol"][:], molsPerRow=len(rOptions["Mol"]), legends=list(rOptions["tag"]))# , subImgSize=(400, 400))
        img.save('Images/test'+str(GroupSelect)+'image.png')
        return rOptions

    def getSets(self, decomp, GroupList, GroupSelect="R1"):
        """function that when given a list of rgoup names will return a pandas dataframe with the corrosponding rgorup data from the file containing data on the rgroups.
        """
        rOptions = pd.DataFrame(columns=["tag", "Smiles"])  # setup dataframe
        if GroupSelect == "R1":  # if reading first R group
            colName = "atag"
            groupnum = 0
        elif GroupSelect == "R2":  # if reading second R group
            colName = "btag"
            groupnum = 1

        for option in range(len(GroupList[groupnum])):
            commonOtherR = decomp.loc[decomp[colName] == GroupList[groupnum][option]]
            commonOtherR = commonOtherR.iloc[0]
            rAtributes = [commonOtherR[colName], commonOtherR[GroupSelect]]
            rOptions.loc[option] = rAtributes
        rdkit.Chem.PandasTools.AddMoleculeColumnToFrame(rOptions, smilesCol='Smiles', molCol='Mol')
        return rOptions

    def makeChoice(self, options, groupname):
        oppdata = pd.DataFrame(options)
        validChoice = False
        while validChoice is False:
            Rchoice = input("> Choice for "+str(groupname)+" " + str(list(oppdata.loc[:,'tag'])) + ": ")                          #add clearer options
            if ("Exit" in Rchoice) or ("exit" in Rchoice) or (Rchoice == "E") or (Rchoice == "e"):
                self.Exit = True
                return
            else:
                self.currentchoice = Rchoice

            print("- you have selected: " + str(self.currentchoice))
            self.chosenmolecule = MolChoose(self.currentchoice[0], self.currentchoice[1])
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
        for groupnum in range(len(self.currentchoice)):
            checkresults = []
            if self.currentchoice[groupnum] in list(self.currentoptions[groupnum].loc[:, "tag"]):
                checkresults.append(True)
            else:
                checkresults.append(False)
                checkcombo = False
                print("\n\n######################## Oh, that looks wrong ########################\n\""+str(self.currentchoice[groupnum])+"\" was not found in the option set of options"+str(list(self.currentoptions[groupnum].loc[:, "tag"]))+",\n please reselect your R Groups\nto exit the game type \"Exit\"\n######################################################################\n")
        self.goodcombo = checkcombo
        return checkcombo

    def UpdateCurrentImage(self):
        from rdkit.Chem.Draw import rdMolDraw2D
        mol = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]  |$;;;;;;;;;;;;R2;;;R1$|')
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)  # or MolDraw2DSVG to get SVGs
        d.drawOptions().addStereoAnnotation = True
        d.DrawMolecule(mol)
        d.FinishDrawing()
        d.WriteDrawingText('Images/CurrentImage.png')

    def chooseAssay(self):
        """
        The user chooses which assay they want to run (currently only pIC50 available).
        The game returns the cost of that assay and lets the user know
        if they have enough money to continue.
        """
        assay_choice = input("\n> Choose assay: ")
        self.assay_cost = assays[assay_choice]['cost']
        print('- The cost of that assay is $' + str(self.assay_cost))
        if self.assay_cost > self.current_balance:
            print('\nYou do not have enough money for another assay')
            self.Exit = True
            return self.endgame()
        else:
            self.assay = assay_choice
            return self.giveFeedback()

    def giveFeedback(self):
        if self.assay == "pIC50":
            pIC50 = self.chosenmolecule.loc[:, "pIC50"].to_string(index=False)
            print("- that molecule's pIC50 is "+str(pIC50))
            self.__scores.append(pIC50)
            # currentstring = self.chosenmolecule.iloc[:,1]
            self.UpdateCurrentImage()
            self.giveNewBalance()

    def endgame(self):
        print("Thank you for playing")
        print(self.__scores)
        self.plotscores()

    def giveStartBalance(self):
        """
        Function returns the starting balance for round 1.
        """
        if self.__loopnumber == 1:
            print("\n Your starting balance is $" + str(self.start_balance) + "\n")
        else:
            pass
    
    def giveNewBalance(self):
        """
        The current balance is updated according to the cost of the assay.
        If the balance reaches 0, the game ends.
        """
        if self.__scores[-1] == ' Not Assayed' or self.__scores[-1] == ' Not Made':
            print("- Assay and/or synthesis was unsuccessful. Your balance is $" +str(self.current_balance))
        else:
            self.current_balance -= self.assay_cost
            print("- Your new balance is $" +str(self.current_balance))
            if self.current_balance == 0:
                print('You have run out of money')
                self.Exit = True
                return self.endgame()
            else:
                pass

    def plotscores(self):
        import numpy as np
        import matplotlib.pyplot as plt
        numericscorelist = []
        fig, ax = plt.subplots()
        for scorenum in range(len(self.__scores)):
            currentscore = self.__scores[scorenum]
            if currentscore == " Not Made":
                colour = "grey"
                currentscore = 0.0
                labelText = "Not Made"
            elif currentscore == " Assay Failed":
                colour = "red"
                currentscore = 0.0
                labelText = "Assay Failed"
            elif currentscore == " Not Assayed":
                colour = "red"
                currentscore = 0.0
                labelText = "Not Assayed"
            elif currentscore == " Inactive":
                colour = "blue"
                currentscore = 0.0
                labelText = "Inactive"

            else:
                currentscore = float(str(currentscore).replace(" ", ''))
                colour = "limegreen"
                labelText = "pIC50"
            ax.scatter(float(scorenum), float(currentscore), c=colour, label=labelText)
            ax.annotate(str(self.choiceHistory[scorenum]).replace("\'", ''), (scorenum, currentscore))
            numericscorelist.append(float(currentscore))

        # ax.legend(False)
        ax.set_xlim(0, len(numericscorelist))
        ax.set_ylim(0, max(numericscorelist))

        plt.title('Chosen Molecule pIC50')
        plt.xlabel("Round")
        plt.ylabel("pIC50")
        plt.show()
        plt.savefig('Images/CurrentImage.png')

    def democheat(self):
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')


game()
