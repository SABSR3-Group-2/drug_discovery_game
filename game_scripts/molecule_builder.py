"""
Make molecules from a scaffold and r groups.

To do:
need to generate self.frags
in update_lead() use ord(self.picked.tag[0]) to join on the [*1] thing
"""
import rdkit
import arcade
import re
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from r_groups_selection import get_selection
from descriptors import get_descriptors

# from inventory import inventory_main

# Cleanse the Images generated in previous rounds
for f_name in os.listdir(os.path.join('Images', 'game_loop_images')):
    if f_name[-4:] == '.png':
        os.remove(os.path.join('Images', 'game_loop_images', f_name))

# Screen title and size
SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
INVENTORY_WIDTH = int(SCREEN_WIDTH / 2.2)
SCREEN_TITLE = "Go Forth and Discover Drugs"

# Constants for sprite scaling
CHARACTER_SCALING = 1
TILE_SCALING = 0.5
MOL_SCALING = 0.5
FILTER_SCALING = 0.3
CURSOR_SCALING = 1


# Calculate all descriptors and store in Dataframe
# data = pd.read_csv('data/r_group_decomp.csv')  # read data
# cols = [c for c in data.columns if re.match('R*\d', c)]  # get the r group cols
# desc = [get_descriptors(x) for x in data[cols[0]].unique()]  # calculate descriptors for all moieties of a given r group
# desc_df = pd.DataFrame(desc)  # make dataframe
# desc_df.insert(0, 'atag', data['atag'].unique())  # insert tags

# Initial scaffold molecule
# scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1] |$;;;;;;;;;;;;R2;;;R1$|')


class MyGame(arcade.Window):
    """
    Main game class
    """

    def __init__(self):

        # Call the parent class and set up the window
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)

        # Initial scaffold molecule
        self.scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')  # |$;;;;;;;;;;;;R2;;;R1$|')

        # Lists that keep track of the 'sprites' aka molecules and r groups
        self.scaffold_list = None
        self.rgroup_list = None
        self.held_molecule = []
        self.scaffold_sprite = None
        arcade.set_background_color(arcade.color.WHITE)

        # Variable to keep track round number
        self.round_count = 0

        """code from inventory.py"""
        # Make relative units for responsive design
        self.vw = int(INVENTORY_WIDTH / 6)  # relative width
        self.vh = int(SCREEN_HEIGHT / 6)  # relative height

        # list to hold the sprites
        self.r_sprite_list = None  # r group sprites
        self.filter_sprite_list = None  # feature filter button sprites
        self.picked_r = None  # the currently held r sprite
        self.current_rs = [0, 0]  # holds the current R1 and R2 groups e.g. [<sprite.object>, <sprite.object>]
        self.lead = self.scaffold  # Initialise with the starting scaffold
        self.tag = 'atag'

        # Track which feature the r groups are being sorted by
        self.feature = None
        self.filters = ['MW', 'logP', 'TPSA', 'HA', 'h_acc', 'h_don', 'rings']

        # Used to keep track of our scrolling
        self.view_top = SCREEN_HEIGHT

        # self.scrolled = 0
        self.top_bound = 0  # the maximum y value
        self.bottom_bound = 0  # the minimum y value

    def _build_lead(self, cur, new, no):
        """
        Helper function for adding a new r group to the current lead:

        :param cur: current lead
        :type cur: `RDkit.mol`
        :param new: the new r group to be added
        :type new: str (SMILES)
        :param no: the number associated with the R group e.g. 'atag' = 1, 'btag' = 2
        :type no: int
        :return: the updated lead
        :rtype: `RDkit.mol`
        """

        mol = Chem.MolToSmiles(cur) + '.' + new
        mol = mol.replace(f'[*:{no}]', '9')
        mol = mol.replace('(9)', '9')
        return Chem.MolFromSmiles(mol)


    def update_lead(self):
        """
        Update the lead with the newly added r groups:
            1. If current = scaffold:
                -> add R group at position
            2. Add R1 to 1*:
                if 1* in lead:
                    -> add R1 to 1*
                else:
                    -> add R1 to .

        c1 = current r1
        c2 = current r2
        ss = starring scaff
        cs = current scaff

        cs = ss + c1 + c2
        cs = ss + c1b + c2

        O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1] = Scaff
        Cc1ccc(S(=O)(=O)NC(=N)NCCC[*:1])cc1 = R1
        O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)9.Cc1ccc(S(=O)(=O)NC(=N)NCCC9)cc1 = SCR1
        Fc1ccc([*:1])cc1 = R1b
        Cc1ccc(S(=O)(=O)NC(=N)NCCCC(NS(=O)(=O)c2ccc([*:2])cc2)C(=O)O)cc1.Fc1ccc9cc1 = SCR1b

        IN: self.picked_r, self.current_rs, self.scaffold
        OUT: self.lead, self.current_rs
        scaffold = self.scaffold
        lead = scaffold
        lead += picked_r
        lead += c for c in current_rs
        self.lead = lead
        """
        r_index = ord(self.picked_r.tag[0].lower()) - 97  # calculate the index given the rtag
        self.current_rs[r_index] = self.picked_r  # update the current_rs with the new pick
        # current_smiles = list(map(lambda picked: self.desc_df.loc[self.desc_df[picked.tag[0]] == picked.tag,
        #                                                           'mol'].item(), self.current_rs))  # sprites -> smiles

        current_scaff = self.scaffold  # initialise with original  scaffold
        for i, c in enumerate(self.current_rs):
            if c != 0:
                current_scaff = self._build_lead(current_scaff, c.smiles, i + 1)
            else:
                pass
        self.lead = current_scaff  # update the current lead

    def make_coordinates(self, n_sprites):
        """Function to make the coordinates for the r sprites.

        :param n_sprites: number of sprites to make coordinates for
        :type n_sprites: int

        :return: list containing 2D coordinates for each sprite, starting in TOP LEFT
        :rtype: nested list
        """

        coordinate_list = []
        full_rows = int(n_sprites / 3)  # number of full rows
        last_row = n_sprites % 3  # number in last, incomplete row

        # make full rows
        for i in range(0, full_rows):
            vh = SCREEN_HEIGHT - (self.vh * (i + 1))  # the y coordinate for a single row with 3 columns from the top
            single_row = [[self.vw, vh],
                          [self.vw * 3, vh],
                          [self.vw * 5, vh]]
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

    def setup_sprites(self, tag, feat='MW'):
        """
        Reads in the specified sprites in a given order e.g. by h_don. Then assigns the coordinates and returns the
        SpriteList object

        :param feat: the feature to sort by
        :type feat: string
        :param tag: the tag of the r group e.g. 'A'
        :type tag: string
        :return: `SpriteList` object
        """
        self.tag = tag
        self.r_sprite_list = arcade.SpriteList(use_spatial_hash=True)
        r_col_index = ord(tag[0].lower()) - 97  # e.g. if 'atag' then 0, 'btag' then 1 etc..

        # Sort r_sprites by specified feature
        data = pd.read_csv('data/r_group_decomp.csv')  # read data
        cols = [c for c in data.columns if re.match('R*\d', c)]  # get the r group cols e.g. ['R1', 'R2']
        desc = [get_descriptors(x) for x in
                data[cols[r_col_index]].unique()]  # calculate descriptors for all moieties of a given r group
        desc_df = pd.DataFrame(desc)  # make dataframe
        desc_df.insert(0, tag, data[tag].unique())  # insert tags
        desc_df.sort_values(feat, inplace=True, ascending=False)  # sort r groups by specified feature
        self.desc_df = desc_df  # for use in key press()
        print(self.desc_df.columns)

        for file in desc_df[tag].unique():
            r_sprite = arcade.Sprite(f'Images/r_group_pngs/{file}.png', MOL_SCALING)
            r_sprite.tag = file  # assign the tag to the r sprite object for later ID
            self.r_sprite_list.append(r_sprite)

        # Create coordinates for the r_sprites
        coordinate_list = self.make_coordinates(len(self.r_sprite_list))

        # Assign coordinates to r_sprites
        for i, sprite in enumerate(self.r_sprite_list):
            sprite.position = coordinate_list[i]

    def setup(self):
        """
        This function sets up the game, call it to restart.
        """
        # update round number
        self.round_count += 1

        # Draw the scaffold molecule
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(self.lead)
        d.FinishDrawing()
        d.WriteDrawingText('Images/game_loop_images/scaffold{}.png'.format(self.round_count))

        # fetch n R groups from the dataset
        # Turn decomp R group dataset into a pandas dataframe
        rgroupdata = pd.read_csv('data/r_group_decomp.csv')

        # Get R1 groups
        self.r1 = get_selection('R1', 3, rgroupdata)

        # get R2 groups
        self.r2 = get_selection('R2', 3, rgroupdata)

        # Create the Rgroup dictionary
        self.rgroup_dict = {}

        # Convert R group SMILES to displayable images
        for r in [self.r1, self.r2]:
            for idx, list in enumerate(r[1]):
                rgroup = Chem.MolFromSmiles(list)
                d = rdMolDraw2D.MolDraw2DCairo(250, 200)
                d.drawOptions().addStereoAnnotation = True
                d.drawOptions().clearBackground = False
                d.DrawMolecule(rgroup)
                d.FinishDrawing()
                d.WriteDrawingText('Images/game_loop_images/rgroup {}.png'.format(r[0][idx]))

        # Create the sprite lists

        self.scaffold_list = arcade.SpriteList()
        self.rgroup_list = arcade.SpriteList(use_spatial_hash=True)

        # Set up the scaffold, placing it at the centre of the screen

        self.scaffold_sprite = arcade.Sprite('Images/game_loop_images/scaffold{}.png'.format(self.round_count),
                                             CHARACTER_SCALING)
        self.scaffold_sprite.position = (int(SCREEN_WIDTH * 0.75), SCREEN_HEIGHT * 0.5)
        self.scaffold_list.append(self.scaffold_sprite)

        # Create and display R2 groups on the left
        # for r in range(len(self.r2[0])):
        #     rgroup = arcade.Sprite('Images/game_loop_images/rgroup {}.png'.format(self.r2[0][r]), TILE_SCALING)
        #     # rgroup.position = (SCREEN_WIDTH / 2, r * 200 + 100)
        #     rgroup.position = (self.vw * 8 , r * 200 + 100)
        #     self.rgroup_list.append(rgroup)
        #     self.rgroup_dict.update({rgroup: f"{self.r2[0][r]}"})

        # Create and display R1 groups on the right
        # for r in range(len(self.r1[0])):
        #     rgroup = arcade.Sprite('Images/game_loop_images/rgroup {}.png'.format(self.r1[0][r]), TILE_SCALING)
        #     rgroup.position = (SCREEN_WIDTH - self.vw * 2, r * 200 + 100)
        #     self.rgroup_list.append(rgroup)
        #     self.rgroup_dict.update({rgroup: f"{self.r1[0][r]}"})

        """
        Code from inventory.py:
            Process to setting up the r_sprites
                1. Create SpriteList of Sprites
                2. Create list of coordinates
                3. Assigning coordinates to each sprite
        """
        # Set up r group sprites
        # Create a list of the sprites
        self.r_sprite_list = arcade.SpriteList(use_spatial_hash=True)

        # Read in the sprite .pngs and create sprites
        self.setup_sprites(self.tag, feat='MW')  # setup sprites with 'atag' as default

        # Set up feature filter buttons
        # Read in filter sprites
        self.filter_sprite_list = arcade.SpriteList(use_spatial_hash=False)
        for i, f in enumerate(self.filters):  # create and position a button for each filter
            filter_sprite = arcade.Sprite(f'Images/filter_pngs/{f}.png', FILTER_SCALING)
            filter_sprite.position = (self.vw * (0.7 * i + 1), SCREEN_HEIGHT - 30)
            filter_sprite.filter = f
            self.filter_sprite_list.append(filter_sprite)

    def on_draw(self):
        """renders the screen"""

        arcade.start_render()

        # Draw sprites
        self.scaffold_list.draw()

        """Code from inventory.py"""
        # Draw r_sprites
        self.r_sprite_list.draw()

        # Draw the menu bar
        arcade.draw_rectangle_filled(INVENTORY_WIDTH / 2,
                                     SCREEN_HEIGHT,
                                     INVENTORY_WIDTH,
                                     self.vh,
                                     color=arcade.color.OXFORD_BLUE)

        arcade.draw_text('press "c" to confirm selection', 0, SCREEN_HEIGHT - int(self.vh / 2) - 20,
                         color=arcade.color.OXFORD_BLUE)

        # Delineate inventory and dragndrop
        arcade.draw_line(INVENTORY_WIDTH, SCREEN_HEIGHT, INVENTORY_WIDTH, 0, arcade.color.OXFORD_BLUE)

        # Draw the filters
        self.filter_sprite_list.draw()

    def on_mouse_press(self, x, y, button, key_modifiers):
        """ Called when the user presses a mouse button. """

        """Code from inventory.py"""
        # Find what the user has clicked on (y coordinate is calculated dynamically)
        clicked = arcade.get_sprites_at_point((x, y), self.filter_sprite_list)

        if len(clicked) > 0:
            # Change the clicked button
            [f._set_color(arcade.color.WHITE) for f in self.filter_sprite_list]
            feature = clicked[0]
            feature._set_color(arcade.color.DARK_CANDY_APPLE_RED)  # turn filter red

            # Sort the r_sprites
            self.desc_df.sort_values(feature.filter, inplace=True, ascending=False)

            # redraw reordered sprites
            self.setup_sprites(self.tag, str(feature.filter))
            self.r_sprite_list.draw()

        # Pick the r sprite and shade it
        for r in self.r_sprite_list:
            r._set_alpha(255)  # removes the shade from any other r_sprites

        self.picked_r = arcade.get_sprites_at_point((x, y), self.r_sprite_list)[-1]  # pick the top sprite
        self.picked_r.smiles = self.desc_df.loc[self.desc_df[self.tag] == self.picked_r.tag, 'mol'].item()  # give smile
        self.picked_r._set_alpha(50)  # shade

    def on_mouse_scroll(self, x: int, y: int, scroll_x: int, scroll_y: int):
        """Redraw the sprites lower instead of scrollling"""

        for i, r in enumerate(self.r_sprite_list):
            r.position = (r.position[0], r.position[1] + scroll_y)
        self.view_top += scroll_y

    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """
        desc_df = self.desc_df
        tag = self.tag
        if symbol == arcade.key.R:
            # Restart
            self.setup()

        if self.picked_r is not None:
            if symbol == arcade.key.C:
                self.update_lead()
                self.setup()
                self.on_draw()

            # if symbol == arcade.key.C:
            #     self.scaffold_list.update()
            #     smiles = desc_df.loc[desc_df[tag] == self.picked_r.tag, 'mol'].item()  # get the smile for the r
            #     # Empty the scaffold lists
            #     print(smiles)
            #     self.scaffold_list = None
            #     self.scaffold_sprite = None
            #     # Combine the initial scaffold and the selected R group into one string
            #     # smiles = smiles + '.' + Chem.MolToSmiles(self.scaffold)
            #     smiles = Chem.MolToSmiles(self.scaffold) + '.' + smiles
            #     # Replace the attachment vector with an integer, if vector is attached to sp2 carbon
            #     new = smiles.replace('([*:1])', '9')
            #     # Replace the attachment vector with an integer, if vector is attached to sp2 carbon
            #     new = new.replace('[*:1]', '9')
            #     print(new)
            #     # Replace scaffold with new combined scaffold. The integers in the smiles indicate the atoms are bonded
            #     self.scaffold = Chem.MolFromSmiles(new)
            #     self.setup()
            #     self.on_draw()
            #
            # if symbol == arcade.key.D:
            #     self.scaffold_list.update()
            #     smiles = desc_df.loc[desc_df[tag] == self.picked_r.tag, 'mol'].item()  # get the smile for the r
            #     # Empty the scaffold lists
            #     print(smiles)
            #     self.scaffold_list = None
            #     self.scaffold_sprite = None
            #     # Combine the initial scaffold and the selected R group into one string
            #     # smiles = smiles + '.' + Chem.MolToSmiles(self.scaffold)
            #     smiles = Chem.MolToSmiles(self.scaffold) + '.' + smiles
            #     # Replace the attachment vector with an integer, if vector is attached to sp2 carbon
            #     new = smiles.replace('([*:2])', '9')
            #     # Replace the attachment vector with an integer, if vector is attached to sp2 carbon
            #     new = new.replace('[*:2]', '9')
            #     print(new)
            #     # Replace scaffold with new combined scaffold. The integers in the smiles indicate the atoms are bonded
            #     self.scaffold = Chem.MolFromSmiles(new)
            #     self.setup()
            #     self.on_draw()

        if symbol == arcade.key.A:
            self.setup_sprites(tag='atag', feat='MW')

        if symbol == arcade.key.B:
            self.setup_sprites(tag='btag', feat='MW')


# Run the game loop
def main():
    """ Main game method"""
    _window = 'Main'
    window = MyGame()
    window.setup()
    arcade.run()


if __name__ == "__main__":
    main()
