"""
Analysis view
"""

import arcade
import pandas as pd
from combine import MolChoose
from rdkit import Chem
from rdkit.Chem import Draw
from math import isnan
from matplotlib import pyplot as plt
import os
from os import listdir
from os.path import isfile, join
from end_game_screen import EndView



SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
MENU_WIDTH = int(SCREEN_WIDTH / 3)
SCREEN_TITLE = "Analysis"

CARD_WIDTH = 130
CARD_HEIGHT = 70
MAT_WIDTH = 150
MAT_HEIGHT = 200

class Card():
    """
    Card class to store molecule info and text coordinates
    """
    def __init__(self, card_coordinates, atag, btag, pic50=None, cl_mouse=None, cl_human=None, logd=None, pampa=None):
        cleaned_input = []
        for assay_value in [pic50, cl_mouse, cl_human, logd, pampa]:
            if type(assay_value) == str:
                cleaned_input.append(assay_value)
            elif isnan(float(assay_value)):
                cleaned_input.append("Not Tested")
            else:
                cleaned_input.append(assay_value)
        self.pic50, self.cl_mouse, self.cl_human, self.logd, self.pampa = cleaned_input


        self.atag = atag
        self.btag = btag
        self.x_tag = card_coordinates[0] - 30
        self.y_tag = card_coordinates[1] + 35
        self.tag_text = f"{self.atag}, {self.btag}"

        self.x_pic50 = card_coordinates[0] - 70
        self.y_pic50 = card_coordinates[1] - 60
        self.pic50_text = f"pIC50: {self.pic50}"

        self.x_cl_mouse = card_coordinates[0] - 70
        self.y_cl_mouse = card_coordinates[1] - 80
        self.cl_mouse_text = f"Cl (mouse): {self.cl_mouse}"

        self.x_cl_human = card_coordinates[0] - 70
        self.y_cl_human = card_coordinates[1] - 100
        self.cl_human_text = f"Cl (human): {self.cl_human}"

        self.x_logd = card_coordinates[0] - 70
        self.y_logd = card_coordinates[1] - 120
        self.logd_text = f"LogD: {self.logd}"

        self.x_pampa = card_coordinates[0] - 70
        self.y_pampa = card_coordinates[1] - 140
        self.pampa_text = f"PAMPA: {self.pampa}"

class ReviewGraph():
    """
    creates graphs for user to review their performace
    """

    def __init__(self, Datatoplot):
        self.data = Datatoplot
        self.fig = plt.figure(figsize=[5, 5])
        self.ax = self.fig.add_subplot(111)
        #self.yvariable = 'pic50'
        #self.xvariable = 'logP'
        self.axisnames = ['','','tags'] 
        self.title = "Review graph"


    def bar(self, yVarr=None, xVar=None, yVar=None):
        pass


    def scatter(self, axisnames, yVarr=None, xVar=None, yVar=None):
        self.axisnames = axisnames
        datavalues = [[],[],[],[]] # [[x],[y],[lables],[colours]] 
        for axisno in [0,1,2]:
            currentaxis =  self.axisnames[axisno]
            if currentaxis in ["logP", "pic50", "cl_mouse", "cl_human", "logd", "pampa", "MW", "logP", "TPSA", "HA", "h_acc", "h_don", "rings"]:
                indexno = 0
                for element in self.data[currentaxis]:
                    if axisno == 0:
                        datavalues[3].append("blue")
                    isfloat = True
                    try:
                        float(element)
                    except ValueError:
                        isfloat = False
                    if isfloat is True:
                        datavalues[axisno].append(float(element))
                    else:
                        datavalues[axisno].append(0.0)
                        if element == "Inactive":
                            datavalues[3][indexno] = ("#909496")
                        elif element == "Not Made":
                            datavalues[3][indexno] = ("#909496")
                        elif element == "Not Assayed":
                            datavalues[3][indexno] = ("#909496")
                        elif element == "NaN":
                            datavalues[3][indexno] = ("#909496")
                        else:
                            datavalues[3][indexno] = ("#909496")
                    indexno += 1
            elif currentaxis in ["tags"]:
                for mol in self.data.iterrows():
                    datavalues[axisno].append(str(mol[1][0])+','+str(mol[1][1]))
        plt.scatter(datavalues[0], datavalues[1], c=datavalues[3])
        for i, datavalues[2] in enumerate(datavalues[2]):
            plt.annotate(datavalues[2], (datavalues[0][i], datavalues[1][i]))
        self.title=str(self.axisnames[1])+" against "+str(self.axisnames[0])
        self.formatgraph()
        plt.savefig(os.path.join("Images","review","maingraph.png"), facecolor='white', transparent=False)
        
        tempdir = os.path.join('Images', 'temp')
        CurFiles = [f for f in listdir(tempdir) if (isfile(join(tempdir, f)) and f.startswith("TempGraph") and f.endswith(".png"))]
        print(f"all: {CurFiles}")
        CurFiles = [s.strip("TempGraph") for s in CurFiles]
        CurFiles = [s.strip(".png") for s in CurFiles]
        print(f"striped: {CurFiles}")
        CurFiles = [int(s) for s in CurFiles if s.isnumeric()]
        print(f"numeric: {CurFiles}")
        if CurFiles != []:
            curMax = max(CurFiles)
        else:
            curMax = 0
        plt.savefig(os.path.join("Images","temp",f"TempGraph{curMax+1}.png"), facecolor='white', transparent=False)
        #return(os.path.join("Images","review",f"maingraph.png"))
        return(os.path.join("Images","temp",f"TempGraph{curMax+1}.png"))


    def formatgraph(self):

        self.ax.set_xlabel(str(self.axisnames[0]))
        self.ax.set_ylabel(str(self.axisnames[1]))
        self.ax.set_title(self.title)



class AnalysisView(arcade.View):
    """
    Analysis view class
    """

    def __init__(self, feedback_view=None):
        super().__init__()
        self.feedback_view = feedback_view
        self.final_df = feedback_view.final_df
        self.button_list = None

        # stores the components of the 'cards'
        self.mol_list = None  # stores the mol sprites (the molecule images on the cards)
        self.mat_list = None  # stores the 'mats' (rectangle representing outside of the card)
        self.text_list = []  # stores the text information for each card (generated with the Card class)

        arcade.set_background_color(arcade.color.WHITE)

        ### GRAPH Init ###
        # stores graphs
        self.graph_list = None
        # store graph related buttons
        self.axisbutton_list = None
        self.axisToggleButton_list = None

        # basic graph properties
        self.buttonscale = 0.4
        self.currentx = "pic50"
        self.currenty = "logP"
        self.axisselectmode = "x"
        self.properties = ["logP", "pic50", "cl_mouse", "cl_human", "logd", "pampa", "MW", "logP", "TPSA", "HA", "h_acc", "h_don", "rings"]

        
        ### Cards Init ##

        # Make relative units for responsive design
        self.vw = int(MENU_WIDTH / 4)  # relative width
        self.vh = int(SCREEN_HEIGHT / 3)  # relative height

        # Used to keep track of our scrolling
        self.view_top = SCREEN_HEIGHT

        # self.scrolled = 0
        self.top_bound = 0  # the maximum y value
        self.bottom_bound = 0  # the minimum y value

        # get font from font file
        self.font = os.path.join('fonts', 'arial.ttf')

        # stores which molecule card has been chosen
        self.mol_choice = None

        self.setup()

    def setup(self):
        self.setupCards()
        self.setupGraph()

    def setupCards(self):
        """
        This function sets up the view, call it to restart.
        """
        ### Card Setup ###
        # create the end button
        self.button_list = arcade.SpriteList()
        end_button = arcade.Sprite(f'Images/button_pngs/end_game_blue.png', 0.5)
        end_button.position = SCREEN_WIDTH - 70, 50
        end_button.name = 'end'
        self.button_list.append(end_button)

        # create a button to allow user to explore a chosen molecule in the molecule_buidler
        builder_button = arcade.Sprite(f'Images/button_pngs/mol_builder.png', 0.5)
        builder_button.position = 50, SCREEN_HEIGHT - 70
        builder_button.name = 'builder'
        self.button_list.append(builder_button)

        # create a button to allow the user to run more assays on a chosen molecule in feedback_buttons
        assays_button = arcade.Sprite(f'Images/button_pngs/run_assays.png', 0.5)
        assays_button.position = 150, SCREEN_HEIGHT - 70
        assays_button.name = 'assays'
        self.button_list.append(assays_button)

        # create the molecule sprites for the cards
        self.mol_list = arcade.SpriteList()
        for index, row in self.feedback_view.final_df.iterrows():
            # get the molecules that have been built/assayed from the final_df
            mol_info = MolChoose(row['atag'], row['btag'], DataSource=os.path.join('data', 'r_group_decomp.csv')).reset_index(drop=True)
            mol = Chem.MolFromSmiles(mol_info.at[0, 'mol'])
            Chem.Draw.MolToFile(mol, 'Images/game_loop_images/selected_mol{}.png'.format(index),
                            size=(CARD_WIDTH, CARD_HEIGHT), imageType=None)  # save image of the molecule to file
            card_sprite = arcade.Sprite('Images/game_loop_images/selected_mol{}.png'.format(index), scale=1)  # create sprites of the molecules
            self.mol_list.append(card_sprite)

        # create blank 'mats' to represent the outline of the cards
        # the mats will be the clickable item for each card to allow the user to select molecules
        self.mat_list = arcade.SpriteList()
        for (index, row), i in zip(self.feedback_view.final_df.iterrows(), range(len(self.mol_list))):
            mat_sprite = arcade.SpriteSolidColor(width = MAT_WIDTH, height = MAT_HEIGHT, color=arcade.color.LIGHT_BLUE)
            # add tag attributes
            mat_sprite.atag = row['atag']
            mat_sprite.btag = row['btag']
            self.mat_list.append(mat_sprite)

        # create coordinates for the mat sprites
        mat_coordinate_list = self.make_coordinates(len(self.mat_list))
        # use the same coordinates but with a higher y value for the mol sprites
        mol_coordinate_list = []
        for coords in mat_coordinate_list:
            mol_coords = []
            mol_coords.append(coords[0])
            mol_coords.append(coords[1] + 50)
            mol_coordinate_list.append(mol_coords)

        # set the coordinates of the molecule sprites
        for i, sprite in enumerate(self.mol_list):
            sprite.position = mol_coordinate_list[i]

        for i, sprite in enumerate(self.mat_list):
            sprite.position = mat_coordinate_list[i]

        # use the Card class to create objects that store the molecule info and coordinates
        for (index, row), coord in zip(self.feedback_view.final_df.iterrows(), mol_coordinate_list):
            cardtext = Card(coord, row['atag'], row['btag'], row['pic50'], row['cl_mouse'], row['cl_human'], row['logd'], row['pampa'])
            self.text_list.append(cardtext)

    def setupGraph(self):
        ### Graph Setup ###
        """
        This function sets up the view, call it to restart.
        """
        # # create end button
        # self.end_button_list = arcade.SpriteList()
        # end_button = arcade.Sprite(os.path.join("Images","button_pngs","end_game_blue.png"), 0.5)
        # end_button.position = SCREEN_WIDTH - 50, 30
        # end_button.name = 'end'
        # self.end_button_list.append(end_button)


        #create graph
        self.cleartempgraphs()
        print("debug a: "+str(self.feedback_view.final_df))
        print("debug a: "+str(self.final_df))
        self.working_graph = ReviewGraph(self.feedback_view.final_df)

        self.graph_list = arcade.SpriteList()
        main_graph = arcade.Sprite(os.path.join('Images','review','maingraph.png'))
        main_graph.position = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH/3)-(SCREEN_WIDTH/20), (SCREEN_HEIGHT*0.4)
        main_graph.name = 'maingraph'
        self.graph_list.append(main_graph)

        
        
        buttonson = True
                


        if buttonson == True:
            #create axis buttons
            self.axisbutton_list = arcade.SpriteList()
            buttonheight = 115*self.buttonscale
            buttonwidth = 85*self.buttonscale
            self.mat_list = arcade.SpriteList()

            # setup button collumns
            numcolls = 7
            colproperties = []
            for col in range(numcolls):
                colproperties.append([])

            for j in range(len(self.properties)):
                colproperties[j%numcolls].append(self.properties[j])
            
            maxlen = []
            for col in colproperties:
                maxlen.append(len(col))
            maxlen = max(maxlen)

            for col in range(numcolls):
                i = 0

                if numcolls == 2:
                    if col == 0:
                        colmodifer = -2.25*buttonwidth
                    else:
                        colmodifer = 2.25*buttonwidth
                
                if numcolls == 3:
                    if col == 0:
                        colmodifer = -2*buttonwidth
                    elif col == 1:
                        colmodifer = 0
                    else:
                        colmodifer = 2*buttonwidth

                
                if numcolls == 5:
                    if col == 0:
                        colmodifer = -5*buttonwidth
                    elif col == 1:
                        colmodifer = -2.5*buttonwidth
                    elif col == 2:
                        colmodifer = 0
                    elif col == 3:
                        colmodifer = 2.5*buttonwidth
                    else:
                        colmodifer = 5*buttonwidth

                if numcolls == 6:
                    if col == 0:
                        colmodifer = -5*buttonwidth
                    elif col == 1:
                        colmodifer = -2.5*buttonwidth
                    elif col == 2:
                        colmodifer = 0
                    elif col == 3:
                        colmodifer = 2.5*buttonwidth
                    elif col == 4:
                        colmodifer = 2.5*buttonwidth
                    else:
                        colmodifer = 5*buttonwidth

                if numcolls == 7:
                    if col == 0:
                        colmodifer = -6*buttonwidth
                    elif col == 1:
                        colmodifer = -4*buttonwidth
                    elif col == 2:
                        colmodifer = -2*buttonwidth
                    elif col == 3:
                        colmodifer = 0
                    elif col == 4:
                        colmodifer = 2*buttonwidth
                    elif col == 5:
                        colmodifer = 4*buttonwidth
                    else:
                        colmodifer = 6*buttonwidth


                for pair in colproperties[col]:
                    hcenter = SCREEN_HEIGHT - ((maxlen* buttonheight) - (i*buttonheight) - 0.5*buttonheight)
                    for side in ["x", "y"]:
                        wcenter = (SCREEN_WIDTH-(SCREEN_WIDTH/3)-1.5*buttonwidth) - buttonwidth + colmodifer
                        property_button = axisButton(pair, self.working_graph, self.buttonscale)
                        property_button.position = wcenter, hcenter
                        self.axisbutton_list.append(property_button)
                        mat_sprite = arcade.SpriteSolidColor(width = int(buttonwidth), height = int(buttonheight), color=arcade.color.LIGHT_BLUE)
                        mat_sprite.position = wcenter,hcenter
                        mat_sprite.value = pair
                        mat_sprite.side = side
                        mat_sprite.equivalant = property_button
                        self.mat_list.append(mat_sprite)
                    i += 1

        if buttonson == True:

            if numcolls == 2:
                colmodifer = -2.25*buttonwidthh
            if numcolls == 3:
                colmodifer = -5*buttonwidth
            if numcolls == 7:
                colmodifer == 9*buttonwidth

            self.axisToggleButton_list = arcade.SpriteList()
            self.xymat_list = arcade.SpriteList()
            hcenter = SCREEN_HEIGHT - (1.5*buttonheight) #(maxlen* buttonheight) + (1* buttonheight) - (0.5*buttonheight)
            for side in ['x','y']:
                self.toggle_button = None
                if side == "x":
                    wcenter = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH*0.32) - buttonwidth + colmodifer + 2.25*buttonwidth
                    state = "on"
                elif side == "y":
                    wcenter = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH*0.32) + buttonwidth + colmodifer + 2.25*buttonwidth
                    state = "off"
                self.toggle_button = xyToggleButton(side, state, self.buttonscale)
                self.toggle_button.position = wcenter, hcenter
                self.axisToggleButton_list.append(self.toggle_button)
                togmat_sprite = arcade.SpriteSolidColor(width = int(buttonwidth), height = int(buttonheight), color=arcade.color.LIGHT_BLUE)
                togmat_sprite.position = wcenter,hcenter
                togmat_sprite.side = side
                togmat_sprite.equivalant = self.toggle_button
                self.xymat_list.append(togmat_sprite)

    def cleartempgraphs(self):
        tempdir = os.path.join('Images', 'temp')
        filestorm = [join(tempdir, f) for f in listdir(tempdir) if (isfile(join(tempdir, f)) and f.startswith("TempGraph") and f.endswith(".png"))]
        for filetorm in filestorm:
            if os.path.exists(filetorm):
                os.remove(filetorm)



    def make_coordinates(self, n_sprites):
        """
        Function to plot the coordinates for the card sprites.
        """
        coordinate_list = []

        n_cards = len(self.feedback_view.final_df)

        full_rows = int(n_sprites / 2)
        last_row = n_sprites % 2

        # make full rows
        for i in range(0, full_rows):
            vh = SCREEN_HEIGHT - (self.vh * (i + 1))  # the y coordinate for a single row with 3 columns from the top
            single_row = [[self.vw, vh],
                          [self.vw * 3, vh]]
            for c in single_row:  # append each coordinate at a time
                coordinate_list.append(c)

        # make last row
        if last_row == 1:
            coordinate_list.append([self.vw, SCREEN_HEIGHT - self.vh * (full_rows + 1)])
        else:
            coordinate_list.append([self.vw, SCREEN_HEIGHT - self.vh * (full_rows + 1)])
            coordinate_list.append([self.vw * 3, SCREEN_HEIGHT - self.vh * (full_rows + 1)])

        # save the maximum and minimum y values of the coordinate list
        self.top_bound = coordinate_list[0][1] + self.vh
        self.bottom_bound = coordinate_list[-1][1] + full_rows
        return coordinate_list

    def on_draw(self):
        """
        Renders the screen
        """
        arcade.start_render()

        # draw the sprites needed for the cards
        self.mat_list.draw()
        self.mol_list.draw()

        # draw the text that goes on the cards
        for c in self.text_list:
            arcade.draw_text(c.tag_text,
                        c.x_tag,
                        c.y_tag,
                        color=arcade.color.BLACK,
                        font_size=10,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.pic50_text,
                        c.x_pic50,
                        c.y_pic50,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.cl_mouse_text,
                        c.x_cl_mouse,
                        c.y_cl_mouse,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.cl_human_text,
                        c.x_cl_human,
                        c.y_cl_human,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.logd_text,
                        c.x_logd,
                        c.y_logd,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.pampa_text,
                        c.x_pampa,
                        c.y_pampa,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')

        # draw the menu bar
        arcade.draw_rectangle_filled(SCREEN_WIDTH/2,
                                     SCREEN_HEIGHT,
                                     SCREEN_WIDTH,
                                     self.vh,
                                     color=arcade.color.OXFORD_BLUE)
        
        arcade.draw_line(SCREEN_WIDTH/3,
                         SCREEN_HEIGHT,
                         SCREEN_WIDTH/3,
                         SCREEN_HEIGHT - (self.vh/2),
                         color=arcade.color.WHITE)
        arcade.draw_line(SCREEN_WIDTH/3,
                         SCREEN_HEIGHT - (self.vh/2),
                         SCREEN_WIDTH/3,
                         0,
                         color=arcade.color.OXFORD_BLUE)

        arcade.draw_text('Select a molecule to investigate further',
                        15,
                        SCREEN_HEIGHT - 25,
                        color=arcade.color.WHITE,
                        font_size=10,
                        font_name=self.font,
                        align='center')

        # draw the button sprites
        self.button_list.draw()


        # draw graph
        #         
        #draw axis buttons
        self.axisbutton_list.draw()
        #draw axis toggle buttons
        self.axisToggleButton_list.draw()

        #draw main graph
    
        self.graph_list[-1].draw()


    def on_mouse_press(self, x, y, button, modifiers):
        """
        Called when the user presses a mouse button. Used for determining what happens
        when the user clicks on a button.
        """
        # check if the user has clicked on a card
        if x < SCREEN_WIDTH/3:
            clicked = arcade.get_sprites_at_point((x, y), self.mat_list)
            if len(clicked) > 0:  # checks a button has been clicked
                [b._set_color(arcade.color.WHITE) for b in self.mat_list]
                choice = clicked[0]
                choice._set_color(arcade.color.YELLOW)  # selected buttons are changed to yellow
                self.mol_choice = [choice.atag, choice.btag]  # record the tags of the chosen molecule

        # check if the user has clicked on a button
        clicked = arcade.get_sprites_at_point((x, y), self.button_list)
        if len(clicked) > 0:  # checks a button has been clicked
            # if end button clicked, csv file created and window closed
            if clicked[0].name == 'end':
                self.feedback_view.final_df.to_csv('data/results.csv', index=False)
                end_view = EndView(self.feedback_view.mol_view)  # create end view and pass mol builder view
                self.window.show_view(end_view)

            # if the molecule builder button is clicked, the chosen molecule tags are passed to self.feedback_view.mol_view
            elif clicked[0].name == 'builder':
                if self.mol_choice is not None:
                    sprites = []
                    # get the sprites in the mol builder script that match the tags
                    for sprite in self.feedback_view.mol_view.r_sprite_list:
                        if sprite.tag == self.mol_choice[0]:
                            sprites.append(sprite)
                        if sprite.tag == self.mol_choice[1]:
                            sprites.append(sprite)
                    # update the lead molecule in mol builder with the correct r groups
                    for sprite in sprites:
                        self.feedback_view.mol_view.picked_r = sprite
                        self.feedback_view.mol_view.picked_r.smiles = self.feedback_view.mol_view.desc_df.loc[self.feedback_view.mol_view.desc_df[self.feedback_view.mol_view.tag] == self.feedback_view.mol_view.picked_r.tag, 'mol'].item()  # give smile
                        self.feedback_view.mol_view.picked_r._set_alpha(50)  # shade
                        self.feedback_view.mol_view.update_lead()
                        self.feedback_view.mol_view.setup()
                        self.feedback_view.mol_view.on_draw()

                    # also update the molecule shown in the feedback view
                    for i, t in enumerate(self.mol_choice):
                        self.feedback_view.tags[i] = t
                    self.feedback_view.setup()
                    self.feedback_view.on_draw()

                    # show the mol builder view
                    self.window.show_view(self.feedback_view.mol_view)

            # if the run assays button is clicked, the chosen molecule tags are again passed to both views
            elif clicked[0].name == 'assays':
                if self.mol_choice is not None:
                    sprites = []
                    # get the sprites in the mol builder script that match the tags
                    for sprite in self.feedback_view.mol_view.r_sprite_list:
                        if sprite.tag == self.mol_choice[0]:
                            sprites.append(sprite)
                        if sprite.tag == self.mol_choice[1]:
                            sprites.append(sprite)
                    # update the lead molecule in mol builder with the correct r groups
                    for sprite in sprites:
                        self.feedback_view.mol_view.picked_r = sprite  # pick the top sprite
                        self.feedback_view.mol_view.picked_r.smiles = self.feedback_view.mol_view.desc_df.loc[self.feedback_view.mol_view.desc_df[self.feedback_view.mol_view.tag] == self.feedback_view.mol_view.picked_r.tag, 'mol'].item()  # give smile
                        self.feedback_view.mol_view.picked_r._set_alpha(50)  # shade
                        self.feedback_view.mol_view.update_lead()
                        self.feedback_view.mol_view.setup()
                        self.feedback_view.mol_view.on_draw()
                    
                    # also update the molecule shown in the feedback view
                    for i, t in enumerate(self.mol_choice):
                        self.feedback_view.tags[i] = t
                    self.feedback_view.setup()
                    self.feedback_view.on_draw()

                    # show the feedback view
                    self.window.show_view(self.feedback_view)
                    arcade.set_background_color(arcade.color.OXFORD_BLUE)
   
        #check axis buttons
        clickedxytoggle = arcade.get_sprites_at_point((x, y), self.xymat_list)
        if len(clickedxytoggle) > 0:  # checks a button has been clicked
            #[b._set_color(arcade.color.WHITE) for b in self.xymat_list]
            xychoice = clickedxytoggle[0]
            for option in self.axisToggleButton_list:
                option.toggle()
            xychoice._set_color(arcade.color.YELLOW)  # selected buttons are changed to yellow
            xyside = xychoice.equivalant.side
            xystate = xychoice.equivalant.state
            if xyside == "x":
                if xystate == "on":
                    self.axisselectmode = "x"
                if xystate == "off":
                    self.axisselectmode = "y"
            elif xyside == "y":
                if xystate == "on":
                    self.axisselectmode = "y"
                if xystate == "off":
                    self.axisselectmode = "x"
            

        #check axis buttons
        clicked = arcade.get_sprites_at_point((x, y), self.mat_list)
        if len(clicked) > 0:  # checks a button has been clicked
            [b._set_color(arcade.color.WHITE) for b in self.mat_list]
            choice = clicked[0]
            choice._set_color(arcade.color.YELLOW)  # selected buttons are changed to yellow
            choiceresult = choice.equivalant.chooseproperty()
            if self.axisselectmode == "x":
                self.currentx = choiceresult[1]
            elif self.axisselectmode == "y":
                self.currenty = choiceresult[1]

            self.plot("scatter")


    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """
        if symbol == arcade.key.SPACE:
            self.plot("scatter")

        if symbol == arcade.key.LEFT:
            self.window.show_view(self.feedback_view)
            arcade.set_background_color(arcade.color.OXFORD_BLUE)
    
    def on_mouse_scroll(self, x: int, y: int, scroll_x: int, scroll_y: int):
        """Redraw the sprites lower instead of scrollling"""

        # change the y position of the sprites
        for i, r in enumerate(self.mol_list):
            r.position = (r.position[0], r.position[1] + scroll_y)
        for i, r in enumerate(self.mat_list):
            r.position = (r.position[0], r.position[1] + scroll_y)
        
        # change the y position of the card text
        for c in self.text_list:
            c.y_tag += scroll_y
            c.y_pic50 += scroll_y
            c.y_cl_mouse += scroll_y
            c.y_cl_human += scroll_y
            c.y_logd += scroll_y
            c.y_pampa += scroll_y

        self.view_top += scroll_y
    

    def plot(self, plottype):

        self.working_graph = ReviewGraph(self.feedback_view.final_df)
        print(f">>>>>>>>>>>>>>>>> {self.feedback_view.final_df}")
        self.graph_list = arcade.SpriteList()
        if plottype == "scatter":
            returnpath = self.working_graph.scatter([self.currentx,self.currenty,"tags"])
            main_graph = arcade.Sprite(returnpath)
            main_graph.position = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH/3)-(SCREEN_WIDTH/20), (SCREEN_HEIGHT*0.4)
            main_graph.name = 'main graph'
            self.graph_list.append(main_graph)
            self.graph_list

class axisButton(arcade.Sprite):
    """Sprite axis button class"""

    def __init__(self, option, graph, scale=0.5):
        # hold the button name and image
        self.button = option
        self.side = "x" #remove this line
        self.image_file_name = os.path.join('Images', 'axisbuttons', f'{self.button}.png')
        self.linkedgraph = graph

        # call the parent class
        super().__init__(self.image_file_name, scale)

    def chooseproperty(self):
        """
        :Returns the assay result when an assay button is clicked.
        :return: assay result for the chosen molecule
        :rtype: string
        """
        print("activating "+str(self.button)+" for "+str(self.side)+" axis")
        # retrieves the appropriate column name
        if self.side == "x":
            return(['x',self.button])
        elif self.side == "y":
            return(['y',self.button])
        else:
            return str(['',self.button])

class xyToggleButton(arcade.Sprite):
    def __init__(self, side, state='off', scale=0.5):
        
        # Set up parent class
        super().__init__(scale=scale)
        
        
        self.side = side
        self.usrscale = scale
        self.state = state



        # Load Textures        
        self.image_file_folder = os.path.join('Images', 'axisbuttons', f'{self.side}')

        #Load textures
        self.texture_off = arcade.load_texture(os.path.join(self.image_file_folder, f'{self.side}_off.png') )
        self.texture_on =  arcade.load_texture( os.path.join(self.image_file_folder, f'{self.side}_on.png') )


        if state == "on":
            self.texture = self.texture_on
        elif state == "off":
            self.texture = self.texture_off

    def toggle(self):
        """
        :Returns the assay result when an assay button is clicked.
        :return: assay result for the chosen molecule
        :rtype: string
        """
        if self.state == 'on':
            self.state = 'off'
        elif self.state == 'off':
            self.state = 'on'
        #self.image_file_name = os.path.join('Images', 'axisbuttons', f'{self.side}_{self.state}.png') 
        self.update(self.state)
        #return([self.side, self.state, self.image_file_name])
    
    def update(self, state=None):
        if state == None:
            state = self.state

        if state == "on":
            self.texture = self.texture_on
        elif state == "off":
            self.texture = self.texture_off


def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = AnalysisView()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()

if __name__ == "__main__":
    main()
