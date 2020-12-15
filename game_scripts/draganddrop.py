"""
Prototype drag and drop molecular constructor for the drug discovery game
"""
import rdkit
import arcade
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from r_groups_selection import get_selection

#Cleanse the Images generated in previous rounds
for f_name in os.listdir('Images'):
    if f_name.startswith('rgroup'):
        os.remove(os.path.join('Images', f_name))

# Screen title and size
SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Go Forth and Discover Drugs"

#Constants for sprite scaling
CHARACTER_SCALING = 1
TILE_SCALING = 0.5

#Scaffold molecule
scaffold = Chem.MolFromSmiles('O=S(C1=CC=CC=C1)(NCC(O)=O)=O')
d = rdMolDraw2D.MolDraw2DCairo(250, 200)
d.drawOptions().addStereoAnnotation = True
d.drawOptions().clearBackground = False
d.DrawMolecule(scaffold)
d.FinishDrawing()
d.WriteDrawingText('Images/scaffold.png')

#fetch n R groups from the dataset
#Turn decomp R group dataset into a pandas dataframe
rgroupdata = pd.read_csv('data/r_group_decomp.csv')
#Get R1 groups
r1 = get_selection('R1', 3, rgroupdata)

#get R2 groups
r2 = get_selection('R2', 3, rgroupdata)

#Convert R1 group SMILES to displayable images 
for r in [r1,r2]:
    print(r)
    for idx, list in enumerate(r[1]):
        rgroup = Chem.MolFromSmiles(list)
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(rgroup)
        d.FinishDrawing()
        d.WriteDrawingText('Images/rgroup {}.png'.format(r[0][idx]))
 
class MyGame(arcade.Window):
    """ 
    Main game class
    """

    def __init__(self):

        #Call the parent class and set up the windwo
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)

        #Lists that keep track of the 'sprites' aka molecules and r groups
        self.scaffold_list = None
        self.rgroup_list = None
        
        self.scaffold_sprite = None

        arcade.set_background_color((255, 255, 255))

    def setup(self):
        """
        This function sets up the game, call it to restart.
        """
        #Create the sprite lists
        self.scaffold_list = arcade.SpriteList()
        self.rgroup_list = arcade.SpriteList(use_spatial_hash=True)

        #Set up the scaffold, placing it at the centre of the screen
        self.scaffold_sprite = arcade.Sprite('Images/scaffold.png', CHARACTER_SCALING )
        self.scaffold_sprite.position = (SCREEN_WIDTH /2,SCREEN_HEIGHT /2)
        self.scaffold_list.append(self.scaffold_sprite)

        #Create and display R1 groups on the left
        for r in range(len(r1[0])):
            rgroup = arcade.Sprite('Images/rgroup {}.png'.format(r1[0][r]), TILE_SCALING)
            rgroup.position = (100, r*200 + 100)
            self.rgroup_list.append(rgroup)
        
        #Create and display R2 groups on the right
        for r in range(len(r2[0])):
            rgroup = arcade.Sprite('Images/rgroup {}.png'.format(r2[0][r]), TILE_SCALING)
            rgroup.position = (800, r*200 + 100)
            self.rgroup_list.append(rgroup)

    

    def on_draw(self):
        #renders the screen

        arcade.start_render()
        
        #Draw sprites
        self.scaffold_list.draw()
        self.rgroup_list.draw()

    def on_mouse_press(self, x, y, button, key_modifiers):
        """ Called when the user presses a mouse button. """

        # Get list of cards we've clicked on
        molecule = arcade.get_sprites_at_point((x, y), self.rgroup_list)

        # Have we clicked on a card?
        if len(molecule) > 0:

            # Might be a stack of cards, get the top one
            primary_molecule = molecule[-1]

            # All other cases, grab the face-up card we are clicking on
            self.held_molecule = [primary_molecule]
            # Save the position
            self.held_molecule_original_position = [self.held_molecule[0].position]
            # Put on top in drawing order
            #self.pull_to_top(self.held_molecule[0])

    def on_mouse_release(self, x: float, y: float, button: int,
                         modifiers: int):
        """ Called when the user presses a mouse button. """
        # If we don't have any cards, who cares
        if len(self.held_molecule) == 0:
            return

        # We are no longer holding cards
        self.held_molecule = []

    def on_mouse_motion(self, x: float, y: float, dx: float, dy: float):
        """ User moves mouse """

        # If we are holding cards, move them with the mouse
        for molecule in self.held_molecule:
            molecule.center_x += dx
            molecule.center_y += dy

    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """
        if symbol == arcade.key.R:
            # Restart
            self.setup()



def main():
    """ Main game method"""
    window = MyGame()
    window.setup()
    arcade.run()

if __name__ == "__main__":
    main()