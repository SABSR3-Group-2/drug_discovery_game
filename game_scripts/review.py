"""
Review view
"""

import arcade
import pandas as pd
from matplotlib import pyplot as plt
import os

SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Review"

class ReviewGraph():
    """
    creates graphs for user to review their performace
    """

    def __init__(self, Datatoplot):
        self.data = Datatoplot
        print("debug (ingraph): "+str(self.data))
        self.fig = plt.figure(figsize=[5, 5])
        self.ax = self.fig.add_subplot(111)
        #self.yvariable = 'pic50'
        #self.xvariable = 'logP'
        self.axisnames = ['','','tags'] 
        self.title = "Review graph"


    def bar(self, yVarr=None, xVar=None, yVar=None):
        pass


    def scatter(self, axisnames, yVarr=None, xVar=None, yVar=None):
        print("axisnames: "+str(axisnames))
        self.axisnames = axisnames
        datavalues = [[],[],[],[]] # [[x],[y],[lables],[colours]] 
        for axisno in [0,1,2]:
            currentaxis =  self.axisnames[axisno]
            print("scatter "+str(self.data))
            if currentaxis in ["logP", "pic50", "cl_mouse", "cl_human", "logd", "pampa", "MW", "logP", "TPSA", "HA", "h_acc", "h_don", "rings"]:
                indexno = 0
                print(self.data.columns)
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
                        print(element)
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
        print(datavalues)
        plt.scatter(datavalues[0], datavalues[1], c=datavalues[3])
        for i, datavalues[2] in enumerate(datavalues[2]):
            plt.annotate(datavalues[2], (datavalues[0][i], datavalues[1][i]))
        self.title=str(self.axisnames[1])+" against "+str(self.axisnames[0])
        self.formatgraph()
        plt.savefig("Images/review/maingraph.png", facecolor='white', transparent=False)


    def formatgraph(self):

        self.ax.set_xlabel(str(self.axisnames[0]))
        self.ax.set_ylabel(str(self.axisnames[1]))
        self.ax.set_title(self.title)


class ReviewView(arcade.View):
    """
    Review class
    """

    def __init__(self, feedback_view=None):
        super().__init__()
        self.feedback_view = feedback_view
        self.final_df = feedback_view.final_df
        self.end_button_list = None
        self.graph_list = None
        self.axisbutton_list = None
        self.axisToggleButton_list = None
        self.buttonscale = 0.4
        self.currentx = "pic50"
        self.currenty = "logP"
        self.axisselectmode = "x"
        arcade.set_background_color(arcade.color.OXFORD_BLUE)

        self.properties = ["logP", "pic50", "cl_mouse", "cl_human", "logd", "pampa", "MW", "logP", "TPSA", "HA", "h_acc", "h_don", "rings"]

        self.setup()


    def setup(self):
        """
        This function sets up the view, call it to restart.
        """
        # create end button
        self.end_button_list = arcade.SpriteList()
        end_button = arcade.Sprite(f'Images/button_pngs/end_game_blue.png', 0.5)
        end_button.position = SCREEN_WIDTH - 50, 30
        end_button.name = 'end'
        self.end_button_list.append(end_button)


        #create graph
        print("debug a: "+str(self.feedback_view.final_df))
        print("debug a: "+str(self.final_df))
        self.working_graph = ReviewGraph(self.feedback_view.final_df)
        # self.working_graph.scatter()

        self.graph_list = arcade.SpriteList()
        main_graph = arcade.Sprite(f'Images/review/maingraph.png')
        main_graph.position = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH/3), (SCREEN_HEIGHT-(SCREEN_HEIGHT/3))/2+(SCREEN_HEIGHT/3)
        main_graph.name = 'maingraph'
        self.graph_list.append(main_graph)

        #self.plot("scatter")
        
        
        buttonson = True
                


        if buttonson == True:
            #create axis buttons
            self.axisbutton_list = arcade.SpriteList()
            buttonheight = 115*self.buttonscale
            buttonwidth = 85*self.buttonscale
            self.mat_list = arcade.SpriteList()

            # setup button collumns
            numcolls = 5
            colproperties = []
            for col in range(numcolls):
                colproperties.append([])
            print(colproperties)

            for j in range(len(self.properties)):
                colproperties[j%numcolls].append(self.properties[j])
            print(colproperties)
            
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


                for pair in colproperties[col]:
                    hcenter = (maxlen* buttonheight) - (i*buttonheight) - 0.5*buttonheight
                    print("col"+str(col))
                    for side in ["x", "y"]:
                        wcenter = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH/3) - buttonwidth + colmodifer
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

            self.axisToggleButton_list = arcade.SpriteList()
            self.xymat_list = arcade.SpriteList()
            hcenter = 0.5*buttonheight #(maxlen* buttonheight) + (1* buttonheight) - (0.5*buttonheight)
            for side in ['x','y']:
                if side == "x":
                    wcenter = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH/3) - buttonwidth + colmodifer + 2.25*buttonwidth
                    state = "on"
                elif side == "y":
                    wcenter = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH/3) + buttonwidth + colmodifer + 2.25*buttonwidth
                    state = "off"
                toggle_button = xyToggleButton(side, state, self.buttonscale)
                toggle_button.position = wcenter, hcenter
                self.axisToggleButton_list.append(toggle_button)
                togmat_sprite = arcade.SpriteSolidColor(width = int(buttonwidth), height = int(buttonheight), color=arcade.color.LIGHT_BLUE)
                togmat_sprite.position = wcenter,hcenter
                togmat_sprite.side = side
                togmat_sprite.equivalant = toggle_button
                self.xymat_list.append(togmat_sprite)

        # self.working_graph = ReviewGraph(self.feedback_view.final_df)

        # self.working_graph.scatter()
        # self.graph_list = arcade.SpriteList()
        # main_graph = arcade.Sprite('Images/review/maingraph.png')
        # main_graph.position = SCREEN_WIDTH/2, SCREEN_HEIGHT/2 
        # main_graph.name = 'main graph'
        # self.graph_list.append(main_graph)
        # #self.graph_list.draw()
        # main_graph.draw()

    def on_draw(self):
        """Renders the screen"""
        arcade.start_render()
        
        # draw the end button
        self.end_button_list.draw()

        #draw axis buttons
        self.axisbutton_list.draw()
        #draw axis toggle buttons
        self.axisToggleButton_list.draw()

        #draw main graph
        self.graph_list.draw()
        

    def on_mouse_press(self, x, y, button, modifiers):
        """
        Called when the user presses a mouse button. Used for determining what happens
        when the user clicks on a button.
        """
        #check axis buttons
        clickedxytoggle = arcade.get_sprites_at_point((x, y), self.xymat_list)
        if len(clickedxytoggle) > 0:  # checks a button has been clicked
            [b._set_color(arcade.color.WHITE) for b in self.xymat_list]
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
            print(self.axisselectmode)

            self.axisToggleButton_list.draw()


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


            print(self.currentx)
            print(self.currenty)

            self.plot("scatter")

            #self.working_graph.scatter()


        # identifies what button the user clicks on
        clicked = arcade.get_sprites_at_point((x, y), self.end_button_list)
        if len(clicked) > 0:  # checks a button has been clicked
            if clicked[0].name == 'end':
                # if end button clicked, csv file created and window closed
                self.feedback_view.final_df.to_csv('data/results.csv', index=False)
                self.window.close()


    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """

        if symbol == arcade.key.SPACE:
            self.plot("scatter")

        if symbol == arcade.key.LEFT:
            self.window.show_view(self.feedback_view)
            arcade.set_background_color(arcade.color.OXFORD_BLUE)


    def plot(self, plottype):

        print("debug b: "+str(self.feedback_view.final_df))
        self.working_graph = ReviewGraph(self.feedback_view.final_df)
        self.graph_list = arcade.SpriteList()
        if plottype == "scatter":
            self.working_graph.scatter([self.currentx,self.currenty,"tags"])

            #main_graph = arcade.Sprite('Images/review/maingraph.png')
            main_graph = arcade.Sprite(f'Images/review/maingraph.png')
            main_graph.position = (SCREEN_WIDTH-(SCREEN_WIDTH/3))/2+(SCREEN_WIDTH/3), (SCREEN_HEIGHT-(SCREEN_HEIGHT/3))/2+(SCREEN_HEIGHT/3)
            main_graph.name = 'main graph'
            self.graph_list.append(main_graph)




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
    """Sprite axis button class"""

    def __init__(self, side, state='off', scale=0.5):
        # hold the button name and image
        self.side = side
        self.state = state
        self.image_file_name = os.path.join('Images', 'axisbuttons', f'{self.side}_{self.state}.png')

        # call the parent class
        super().__init__(self.image_file_name, scale)

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
        self.image_file_name = os.path.join('Images', 'axisbuttons', f'{self.side}_{self.state}.png')   
        self.draw()
        return([self.side, self.state])






def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = ReviewView()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()


if __name__ == "__main__":
    main()
