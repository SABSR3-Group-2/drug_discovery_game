"""
Analysis view
"""

import arcade
import pandas as pd
from matplotlib import pyplot as plt

SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Analysis"

class AnalysisGraph():
    """
    creates graphs for user to review their performace
    """

    def __init__(self, Datatoplot):
        self.data = Datatoplot
        self.fig = plt.figure(figsize=[5, 5])
        self.ax = self.fig.add_subplot(111)
        self.yvariable = 'pIC50'
        self.xvariable = 'tags'

    def bar(self, yVarr=None, xVar=None, yVar=None):

        if type(yVar) == str:
            print(self.yvariable)
            print(yVar)
            self.yvariable = yVar
            print(self.yvariable)
        if type(xVar) == str:
            self.xvariable = xVar

        xvalues = []
        if self.xvariable == "tags":
            for mol in self.data.iterrows():
                print(mol)
                print(mol[1][1])
                xvalues.append(str(mol[1][0])+','+str(mol[1][1]))



        yvalues = self.data[self.yvariable]
        #canvas = agg.FigureCanvasAgg(fig)
        # print(xvalues)
        # print(yvalues)
        plt.bar(xvalues,yvalues)
        self.title=str(self.yvariable)+" against "+str(self.xvariable)
        self.formatgraph()
        plt.savefig("Images/analysis/maingraph.png", facecolor='white', transparent=False)

    def formatgraph(self):


        self.ax.set_ylabel(str(self.yvariable))
        self.ax.set_xlabel(str(self.xvariable))
        self.ax.set_title(self.title)





df = pd.DataFrame({'ATag': ["a10", "a11", "a12"], 'BTag': ['b1', 'c', 'k'], 'pIC50': [100, 110, 120]})

graph = AnalysisGraph(df)
graph.bar('pIC50')

class AnalysisView(arcade.View):
    """
    Analysis class
    """

    def __init__(self, feedback_view=None):
        super().__init__()
        self.feedback_view = feedback_view
        self.end_button_list = None

        arcade.set_background_color(arcade.color.WHITE)

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

    def on_draw(self):
        """Renders the screen"""
        arcade.start_render()

        # draw the end button
        self.end_button_list.draw()

    def on_mouse_press(self, x, y, button, modifiers):
        """
        Called when the user presses a mouse button. Used for determining what happens
        when the user clicks on a button.
        """

        # identifies what button the user clicks on
        clicked = arcade.get_sprites_at_point((x, y), self.end_button_list)
        if len(clicked) > 0:  # checks a button has been clicked
            if clicked[0].name == 'end':
                # if end button clicked, csv file created and window closed
                self.feedback_view.final_df.to_csv('data/results.csv', index=False)
                self.window.close()

    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """
        # df containing results printed when user presses SPACE
        if symbol == arcade.key.SPACE:
            print(self.feedback_view.final_df)

        if symbol == arcade.key.LEFT:
            self.window.show_view(self.feedback_view)
            arcade.set_background_color(arcade.color.OXFORD_BLUE)


def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = AnalysisView()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()


if __name__ == "__main__":
    main()
