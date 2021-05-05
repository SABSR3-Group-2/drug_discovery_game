import os
import arcade
from game_scripts.combine import MolChoose
from game_scripts.descriptors import get_descriptors
from game_scripts.filters import compound_check
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import global_vars

"""
Window for ending the game - displays a comparison of chosen molecule vs final target molecule
"""

# Constants
SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "The Final Selection"

CARD_WIDTH = 250
CARD_HEIGHT = 400
MAT_WIDTH = 250
MAT_HEIGHT = 650


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

        # stores the path to the font file
        self.font = os.path.join('fonts', 'arial.ttf')

        # sets the background color
        arcade.set_background_color(arcade.color.WHITE)

        self.setup()

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
        left_mat.position = 350, 400
        self.mat_list.append(left_mat)

        right_mat = arcade.SpriteSolidColor(CARD_WIDTH, CARD_HEIGHT, arcade.color.LIGHT_BLUE)
        right_mat.position = 650, 400
        self.mat_list.append(right_mat)

        #Generate image of the selected final molecule
        self.card_list = arcade.SpriteList()
        mol_info = MolChoose(self.analysis_view.mol_choice[0], self.analysis_view.mol_choice[1], DataSource=os.path.join('data', 'r_group_decomp.csv')).reset_index(drop=True)
        mol = Chem.MolFromSmiles(mol_info.at[0, 'mol'])
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(mol)
        d.FinishDrawing()
        d.WriteDrawingText('Images/game_loop_images/final_mol.png')
        self.chosen_sprite = arcade.Sprite('Images/game_loop_images/final_mol.png')
        self.chosen_sprite.position = (350, 520)
        self.card_list.append(self.chosen_sprite)

    def on_draw(self):
        """Render the screen"""

        # clear the screen to the background colour
        arcade.start_render()

        self.button_list.draw()

        #draw the 'card mats'
        self.mat_list.draw()
        self.card_list.draw()

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
