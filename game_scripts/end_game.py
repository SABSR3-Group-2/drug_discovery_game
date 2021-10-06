import os
import numpy as np
import arcade
from combine import MolChoose
from descriptors import get_descriptors
from filters import compound_check
import textwrap
from feedback_dict import feedback_clearance, feedback_lipophilicity, feedback_pic50
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
CARD_HEIGHT = 310
MAT_WIDTH = 250
MAT_HEIGHT = 650
MAT_Y_COORD = 530
LEFT_MAT_X_COORD = 350
RIGHT_MAT_X_COORD = 650

CORE = 'O=C(O)C(NS(=O)(=O)c1ccccc1)'

#Constants for spider plot
labels=['pIC50', 'logD', 'Mouse Clearance', 'Human Clearance', 'Permeability']
markers = [0, 1, 2, 3, 4, 5, 6, 7, 8]
str_markers = ["0", "1", "2", "3", "4", "5", '6', '7', '8']
target_data = [7.7, 1.08, 1, 1, 7]

#Dictionaries converting non-numerical assay results into numerical approximate values
clearance_dict = {
    'low (< 5.6)': 1,
    'medium (5.6-30.5)': 4,
    'low (< 3.7)': 1, 
    'good':1,
    'high (> 30.5)': 7,
    'fair': 4,
    'poor': 7,
    'low (< 12)': 1,
    'medium (12-44)': 4,
    'medium (5.6-30.5)':4
}
pampa_dict = {
    'neg':0,
    'poor':1,
    'low': 2.5,
    'fair':5.5,
    'med2high':5.5,
    'good':6.5,
    'best':8
}


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
        self.choice_text_list = []
        self.target_text_list = []
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

        fig= plt.figure(figsize=(6,3.5))
        ax = fig.add_subplot(111, polar=True)
        ax.plot(angles, stats, 'o-', linewidth=2, label="Your Choice")
        ax.fill(angles, stats, alpha=0.25)
        ax.plot(angles, target_data, 'o-', linewidth=2, label='Final Compound')
        ax.fill(angles, target_data, alpha=0.25)
        ax.set_thetagrids(angles[:-1] * 180/np.pi, labels)
        # Go through labels and adjust alignment based on where
        # it is in the circle.
        for label, angle in zip(ax.get_xticklabels(), angles):
            if angle == 0:
                label.set_horizontalalignment('left')
            elif angle in (2.5132741228718345, 3.7699111843077517):
                label.set_horizontalalignment('right')
            else:
                label.set_horizontalalignment('center')

        plt.yticks(plot_markers, plot_str_markers)
        ax.legend(loc='best', bbox_to_anchor=(0.4, 0))
        ax.grid(True)
        plt.tight_layout()

        fig.savefig('images/game_loop_images/spider_plot.png', transparent=True)

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
        fig.savefig('images/game_loop_images/pic50_line.png', transparent=True)

    def setup(self):
        """
        This function sets up the view
        """

         # create end button
        self.button_list = arcade.SpriteList()
        end_button = arcade.Sprite(f'images/button_pngs/end_game_blue.png', 0.5)
        end_button.position = SCREEN_WIDTH/2 , 50
        end_button.name = 'end'
        self.button_list.append(end_button)

        #Create the placeholders for the two 'card mats'
        self.mat_list = arcade.SpriteList()

        left_mat = arcade.SpriteSolidColor(CARD_WIDTH, CARD_HEIGHT, arcade.color.OXFORD_BLUE)
        left_mat.position = LEFT_MAT_X_COORD, MAT_Y_COORD
        self.mat_list.append(left_mat)

        left_background = arcade.SpriteSolidColor(CARD_WIDTH - 10, 150, arcade.color.WHITE)
        left_background.position = LEFT_MAT_X_COORD, 545
        self.mat_list.append(left_background)

        right_mat = arcade.SpriteSolidColor(CARD_WIDTH, CARD_HEIGHT, arcade.color.OXFORD_BLUE)
        right_mat.position = RIGHT_MAT_X_COORD, MAT_Y_COORD
        self.mat_list.append(right_mat)

        right_background = arcade.SpriteSolidColor(CARD_WIDTH - 10, 150, arcade.color.WHITE)
        right_background.position = RIGHT_MAT_X_COORD, 545
        self.mat_list.append(right_background)

        #Align Roche choice and user choice
        core = Chem.MolFromSmiles(CORE)
        target = Chem.MolFromSmiles('CCSCC[C@H](NS(=O)(=O)c1ccc(cc1)c2ccc(SC)cc2)C(=O)O')
        mol_info = MolChoose(self.analysis_view.mol_choice[0], self.analysis_view.mol_choice[1], DataSource=os.path.join('data', 'r_group_decomp.csv')).reset_index(drop=True)
        smiles = mol_info.at[0, 'mol']
        mol = Chem.MolFromSmiles(mol_info.at[0, 'mol'])
        AllChem.Compute2DCoords(core)
        for m in [target, mol]:
            _ = AllChem.GenerateDepictionMatching2DStructure(m, core)


        #Generate sprite for the Roche choice
        self.card_list = arcade.SpriteList()
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(target)
        d.FinishDrawing()
        d.WriteDrawingText('images/game_loop_images/target_mol.png')
        self.target_sprite = arcade.Sprite('images/game_loop_images/target_mol.png')
        self.target_sprite.position = (RIGHT_MAT_X_COORD, MAT_Y_COORD + 20)
        self.card_list.append(self.target_sprite)


        #Generate image of the selected final molecule
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(mol)
        d.FinishDrawing()
        d.WriteDrawingText('images/game_loop_images/final_mol.png')
        self.chosen_sprite = arcade.Sprite('images/game_loop_images/final_mol.png')
        self.chosen_sprite.position = (LEFT_MAT_X_COORD, MAT_Y_COORD + 20)
        self.card_list.append(self.chosen_sprite)

        #Get descriptor information for the final chosen molecule
        descriptor_dict = get_descriptors(smiles)
        if mol_info.at[0, 'pic50'] not in ['Assay Failed', 'Inactive', 'Not Made']:
            pic50 = float(mol_info.at[0, 'pic50'])
        else:
            pic50 = 0
        logd = float(mol_info.at[0, 'logd'])
        mouse_clearance = clearance_dict[mol_info.at[0, 'clearance_mouse']]
        human_clearance = clearance_dict[mol_info.at[0, 'clearance_human']]
        permeability = pampa_dict[mol_info.at[0, 'pampa']]
        
        #List of parameters for final molecule spider plot
        final_property_list = [pic50, logd, mouse_clearance, human_clearance, permeability]

        #Generate spider plot
        self.make_radar_chart(stats=final_property_list)
        self.radarplot = arcade.Sprite('images/game_loop_images/spider_plot.png')
        self.radarplot.position = (250, 200)
        self.card_list.append(self.radarplot)

        #generate text list
        pic50 = mol_info.at[0, 'pic50']
        self.choice_text_list = [f'pIC50: {mol_info.at[0, "pic50"]}', f'logD: {logd}', f'Mouse clearance: {mol_info.at[0, "clearance_mouse"]}', f'Human clearance: {mol_info.at[0, "clearance_human"]}', f'Permeability: {mol_info.at[0, "pampa"]}']
        self.target_text_list = ['pIC50: 7.7', 'logD: 1.08', 'Mouse clearance: low (<5.6)', 'Human clearance: low (<12)', 'Permeability: medium - high']
        
        #Generate feedback text
        pic50 = float(pic50)
        if pic50 < 6.5:
            pic50_feedback = 'low'
        else:
            pic50_feedback = 'good'

        if 0.95 < logd <1.15:
            logd_feedback = 'good'
        elif logd < 0.95:
            logd_feedback = 'low'
        else:
            logd_feedback = 'high'

        if human_clearance != 1:
            clearance_feedback = 'high'
        else:
            clearance_feedback = 'low'

        self.feedback = [feedback_pic50[pic50_feedback], feedback_lipophilicity[logd_feedback], feedback_clearance[clearance_feedback]]
        

    def on_draw(self):
        """Render the screen"""

        # clear the screen to the background colour
        arcade.start_render()

        self.button_list.draw()

        #draw the 'card mats'
        self.mat_list.draw()
        self.card_list.draw()

        #Drawing box titles
        arcade.draw_text('Your Choice', 300, 628, color=arcade.color.WHITE, font_size=15)
        arcade.draw_text("Roche's Choice", 590, 628, color=arcade.color.WHITE, font_size=15)

        #Draw in choice Data
        for counter, item in enumerate(self.choice_text_list):
            arcade.draw_text(item, 240, 450 - counter*17, color=arcade.color.WHITE, font_size=10)

        #Draw in target data
        for counter, item in enumerate(self.target_text_list):
            arcade.draw_text(item, 545, 450 - counter*17, color=arcade.color.WHITE, font_size=10)

        #Draw in feedback text
        feedback = ""
        for item in self.feedback:
            string = textwrap.fill(item, 70)
            feedback += string + '\n' + '\n'
        arcade.draw_text(feedback, 500, 100, color=arcade.color.BLACK)

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
