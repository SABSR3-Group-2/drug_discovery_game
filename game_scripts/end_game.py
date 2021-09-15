import os
import numpy as np
import arcade
from game_scripts.combine import MolChoose
from game_scripts.descriptors import get_descriptors
from game_scripts.filters import compound_check
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import global_vars
import matplotlib.pyplot as plt

"""
Window for ending the game - displays a comparison of chosen molecule vs final target molecule
"""

# Constants
SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "The Final Selection"

CARD_WIDTH = 250
CARD_HEIGHT = 200
MAT_WIDTH = 250
MAT_HEIGHT = 650

CORE = 'O=C(O)C(NS(=O)(=O)c1ccccc1)'

#Constants for spider plot
labels=['pIC50', 'logP', 'logD', 'H-Bond Acceptors', 'H-Bond Donors']
markers = [0, 1, 2, 3, 4, 5, 6, 7, 8]
str_markers = ["0", "1", "2", "3", "4", "5", '6', '7', '8']
target_data = [7.7, 3.95, 1.08, 5, 2]

class EndGame(arcade.View):
    """
    Main class for endgame view
    """

    def __init__(self, analysis_view = None):
        #call the parent class and set up the window
        super().__init__()
        self.analysis_view = analysis_view
        self.button_list = None
        self.card_list = None
        self.mat_list = None
        self.chosen_sprite = None
        self.numberline = None
        self.target_sprite = None
        self.text_list = []
        self.radarplot = None

        # stores the path to the font file
        self.font = os.path.join('fonts', 'arial.ttf')

        # sets the background color
        arcade.set_background_color(arcade.color.WHITE)

        self.setup()

    def make_radar_chart(self, stats, attribute_labels = labels, plot_markers = markers, plot_str_markers = str_markers, target_data = target_data):

        labels = np.array(attribute_labels)

        angles = np.linspace(0, 2*np.pi, len(labels), endpoint=False)
        #The below concatenations are needed to ensure that the input lists have the same lengths, in the final product the last entry in the list must be the same as the first
        stats = np.concatenate((stats,[stats[0]]))
        angles = np.concatenate((angles,[angles[0]]))

        
        target_data = np.concatenate((target_data, [target_data[0]]))

        fig= plt.figure(figsize=(6,3))
        ax = fig.add_subplot(111, polar=True)
        ax.plot(angles, stats, 'o-', linewidth=2, label="GSK's Choice")
        ax.fill(angles, stats, alpha=0.25)
        ax.plot(angles, target_data, 'o-', linewidth=2, label='Final Choice')
        ax.fill(angles, target_data, alpha=0.25)
        ax.set_thetagrids(angles[:-1] * 180/np.pi, labels)
        plt.yticks(plot_markers, plot_str_markers)
        ax.legend(loc='best', bbox_to_anchor=(1, 0.4))
        ax.grid(True)

        fig.savefig('Images/game_loop_images/spider_plot.png', transparent=True)

        return


    def make_numberline(self, mol_info):
        #Generate number-line of pIC50's for comparison with final plot
        # set up the figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(0,10)
        ax.set_ylim(0,10)
        ax.set_title('pIC50', y=0.1, pad=-14)

        # draw lines
        xmin = 0
        xmax = 8
        y = 5
        height = 1

        plt.hlines(y, xmin, xmax)
        plt.vlines(xmin, y - height / 2., y + height / 2.)
        plt.vlines(xmax, y - height / 2., y + height / 2.)

        if mol_info.at[0, 'pic50'] != 'Not Made':
            # draw a point on the line
            px = mol_info.at[0, 'pic50']
            plt.plot(px,y, 'ro', ms = 15, mfc = 'r')
        elif mol_info.at[0, 'pic50'] == 'Inactive':
            px = 0
            plt.plot(px, y, 'ro', ms = 15, mfc = 'r')

        # add numbers
        plt.text(xmin - 0.1, y, '0.0', horizontalalignment='right')
        plt.text(xmax + 0.1, y, '8.0', horizontalalignment='left')
        

        plt.axis('off')
        fig.savefig('Images/game_loop_images/pic50_line.png', transparent=True)

    def setup(self):
        """
        This function sets up the view
        """

         # create end button
        self.button_list = arcade.SpriteList()
        end_button = arcade.Sprite(f'Images/button_pngs/end_game_blue.png', 0.5)
        end_button.position = SCREEN_WIDTH/2 , 50
        end_button.name = 'end'
        self.button_list.append(end_button)

        #Create the placeholders for the two 'card mats'
        self.mat_list = arcade.SpriteList()

        left_mat = arcade.SpriteSolidColor(CARD_WIDTH, CARD_HEIGHT, arcade.color.LIGHT_BLUE)
        left_mat.position = 350, 500
        self.mat_list.append(left_mat)

        right_mat = arcade.SpriteSolidColor(CARD_WIDTH, CARD_HEIGHT, arcade.color.LIGHT_BLUE)
        right_mat.position = 650, 500
        self.mat_list.append(right_mat)

        #Align GSK choice and user choice
        core = Chem.MolFromSmiles(CORE)
        target = Chem.MolFromSmiles('CCSCC[C@H](NS(=O)(=O)c1ccc(cc1)c2ccc(SC)cc2)C(=O)O')
        mol_info = MolChoose(self.analysis_view.mol_choice[0], self.analysis_view.mol_choice[1], DataSource=os.path.join('data', 'r_group_decomp.csv')).reset_index(drop=True)
        smiles = mol_info.at[0, 'mol']
        mol = Chem.MolFromSmiles(mol_info.at[0, 'mol'])
        AllChem.Compute2DCoords(core)
        for m in [target, mol]:
            _ = AllChem.GenerateDepictionMatching2DStructure(m, core)


        #Generate sprite for the GSK choice
        self.card_list = arcade.SpriteList()
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(target)
        d.FinishDrawing()
        d.WriteDrawingText('Images/game_loop_images/target_mol.png')
        self.target_sprite = arcade.Sprite('Images/game_loop_images/target_mol.png')
        self.target_sprite.position = (650, 520)
        self.card_list.append(self.target_sprite)


        #Generate image of the selected final molecule
        d = rdMolDraw2D.MolDraw2DCairo(250, 520)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(mol)
        d.FinishDrawing()
        d.WriteDrawingText('Images/game_loop_images/final_mol.png')
        self.chosen_sprite = arcade.Sprite('Images/game_loop_images/final_mol.png')
        self.chosen_sprite.position = (350, 520)
        self.card_list.append(self.chosen_sprite)

        #Get descriptor information for the final chosen molecule
        descriptor_dict = get_descriptors(smiles)
        if mol_info.at[0, 'pic50'] not in ['Assay Failed', 'Inactive', 'Not Made']:
            pic50 = float(mol_info.at[0, 'pic50'])
        else:
            pic50 = 0
        logd = float(mol_info.at[0, 'logd'])
        
        #List of parameters for final molecule spider plot
        final_property_list = [pic50, descriptor_dict['logP'], logd, descriptor_dict['h_acc'], descriptor_dict['h_don']]
        print(final_property_list)

        #Generate spider plot
        self.make_radar_chart(stats=final_property_list)
        self.radarplot = arcade.Sprite('Images/game_loop_images/spider_plot.png')
        self.radarplot.position = (480, 250)
        self.card_list.append(self.radarplot)

        #generate text list
        pic50 = mol_info.at[0, 'pic50']
        self.text_list = ['pIC50: {}'.format(pic50), 'pIC50: 7.7']

    def on_draw(self):
        """Render the screen"""

        # clear the screen to the background colour
        arcade.start_render()

        self.button_list.draw()

        #draw the 'card mats'
        self.mat_list.draw()
        self.card_list.draw()

        #Drawing box titles
        arcade.draw_text('Your Choice', 310, 580, color=arcade.color.BLACK)
        arcade.draw_text("GSK's Choice", 610, 580, color=arcade.color.BLACK)

        #Draw in pIC50 Data
        arcade.draw_text(self.text_list[0], 240, 420, color=arcade.color.BLACK,
                            font_size=20)
        arcade.draw_text(self.text_list[1], 540, 420, color=arcade.color.BLACK, font_size=20)

    def on_mouse_press(self, x, y, button, modifiers):
        """
        Called when the user presses a mouse button. Used for determining what happens
        when the user clicks on a button
        """

        #Checks to see if the End button has been pressed
        clicked = arcade.get_sprites_at_point((x,y), self.button_list)
        if len(clicked) > 0:
            if clicked[0].name == 'end':
                self.window.close()
        
        

def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = EndGame()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()


if __name__ == "__main__":
    main()
