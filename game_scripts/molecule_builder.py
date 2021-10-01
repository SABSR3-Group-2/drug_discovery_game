"""
Make molecules from a scaffold and r groups.

To do:
"""
import arcade
import re
import os
import pandas as pd
import time
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.Draw import rdMolDraw2D
from descriptors import get_descriptors
from feedback_buttons import FeedbackView
import textwrap

# Cleanse the images generated in previous rounds
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


class MolView(arcade.View):
    """
    Main game class
    """

    def __init__(self, feedback_view=None):

        # Call the parent class and set up the window
        super().__init__()
        self.feedback_view = feedback_view
        RDLogger.DisableLog('rdApp.*')  # don't print warning message for hydrogen removal
        # Initial scaffold molecule
        self.scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')  # |$;;;;;;;;;;;;R2;;;R1$|')
        self.num_vecs = self.get_num_vectors()  # get an integer value for how many vectors there are

        # Lists that keep track of the 'sprites' aka molecules and r groups
        self.scaffold_list = None
        self.scaffold_sprite = None
        arcade.set_background_color(arcade.color.WHITE)

        # Variable to keep track round number
        self.round_count = 0

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
        self.buttons = None  # holds the button sprites for changing between different inventories

        # Track which feature the r groups are being sorted by
        self.feature = 'MW'
        self.filters = ['MW', 'logP', 'TPSA', 'HA', 'h_acc', 'h_don', 'rings']

        # Used to keep track of our scrolling
        self.view_top = SCREEN_HEIGHT

        # Track which sprite we're near for displaying help
        self.hovered = None
        self.hover_time = 0
        self.location = (0, 0)
        self.display_hover = False

        # self.scrolled = 0
        self.top_bound = 0  # the maximum y value
        self.bottom_bound = 0  # the minimum y value

        # create df to record assay data
        self.col_names = [
            'atag', 'btag',
            'pic50', 'cl_mouse', 'cl_human', 'logd', 'pampa',
            'MW', 'logP', 'TPSA', 'HA', 'h_acc', 'h_don', 'rings'
        ]
        self.assay_df = pd.DataFrame(columns=self.col_names)

    def get_num_vectors(self):
        """Deduces the number of vectors on the scaffold molecule from the scaffold smile"""
        vecs = re.findall('\[\*\:\d+\]', Chem.MolToSmiles(self.scaffold))
        return len(vecs)

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

        for file in desc_df[tag].unique():
            r_sprite = arcade.Sprite(f'Images/r_group_pngs/{file}.png', MOL_SCALING)
            r_sprite.tag = file  # assign the tag to the r sprite object for later ID
            self.r_sprite_list.append(r_sprite)

        # Create coordinates for the r_sprites
        coordinate_list = self.make_coordinates(len(self.r_sprite_list))

        # Assign coordinates to r_sprites
        for i, sprite in enumerate(self.r_sprite_list):
            sprite.position = coordinate_list[i]

        # self.setup_values()

    def draw_values(self):
        """
        For the r sprite inventory, this function draws the value that the list is being sorted by under each sprite.
        E.g. if sorting by rings and the R group has 3 rings it will print "3" under the associated sprite.
        """

        for i, r in enumerate(self.desc_df[self.feature]):
            coord = self.r_sprite_list[i].position
            value_coord = [f'{float(r):.0f}', coord[0], coord[1] - 55]
            arcade.draw_text(value_coord[0], value_coord[1], value_coord[2],
                             color=arcade.color.BLACK, align="center", font_size=11)

    def draw_hover(self):

        # Specify the help texts
        text_dict = {'MW': 'Molecular Weight',
                     'logP': 'Lipophilicity',
                     'TPSA': 'Total Polar Surface Area',
                     'HA': 'Heavy Atom count',
                     'h_acc': 'hydrogen bond acceptor count',
                     'h_don': 'hydrogen bond donor count',
                     'rings': 'ring count',
                     'left': 'previous R groups',
                     'right': 'next R groups'}

        text = text_dict[self.hovered.tag]  # what text to write out
        loc = (self.hovered.position[0] + 30, self.hovered.position[1] - 30)  # where to draw it

        # Create the text sprite
        text_sprite = arcade.draw_text(text, loc[0], loc[1], color=arcade.color.BLACK, font_size=10)

        # Draw the background
        width = text_sprite.width
        height = text_sprite.height
        arcade.draw_rectangle_filled(loc[0] + width * 0.5, loc[1] + height * 0.5,
                                     width + 10, height + 10,
                                     color=arcade.color.YELLOW)

        # Draw the text
        text_sprite.draw()

    def combine_sprite_lists(self, lists):
        """Helper method for combining multiple sprite lists into one sprite list"""
        merged = arcade.SpriteList(use_spatial_hash=False)
        for l in lists:
            for s in l:
                merged.append(s)
        return merged

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

        # Create the sprite lists
        self.scaffold_list = arcade.SpriteList()

        # Set up the scaffold, placing it at the centre of the screen
        self.scaffold_sprite = arcade.Sprite('Images/game_loop_images/scaffold{}.png'.format(self.round_count),
                                             CHARACTER_SCALING)
        self.scaffold_sprite.position = (int(SCREEN_WIDTH * 0.75), SCREEN_HEIGHT * 0.5)
        self.scaffold_list.append(self.scaffold_sprite)

        # Set up r group sprites
        # Create a list of the sprites
        self.r_sprite_list = arcade.SpriteList(use_spatial_hash=True)

        # Read in the sprite .pngs and create sprites
        self.setup_sprites(self.tag, feat=self.feature)  # setup sprites with 'atag' as default

        # Set up feature filter buttons
        # Read in filter sprites
        self.filter_sprite_list = arcade.SpriteList(use_spatial_hash=False)
        for i, f in enumerate(self.filters):  # create and position a button for each filter
            filter_sprite = arcade.Sprite(f'Images/filter_pngs/{f}.png', FILTER_SCALING)
            filter_sprite.position = (self.vw * (0.7 * i + 1), SCREEN_HEIGHT - 30)
            filter_sprite.tag = f
            self.filter_sprite_list.append(filter_sprite)

        # Set up inventory navigation button sprites
        self.buttons = arcade.SpriteList(use_spatial_hash=True)
        button = arcade.Sprite(os.path.join('Images', 'filter_pngs', 'left_arrow.png'), FILTER_SCALING / 2)
        button.position = (self.vw * 0.3, SCREEN_HEIGHT - 30)
        button.tag = 'left'
        self.buttons.append(button)
        button = arcade.Sprite(os.path.join('Images', 'filter_pngs', 'right_arrow.png'), FILTER_SCALING / 2)
        button.position = (INVENTORY_WIDTH - self.vw * 0.3, SCREEN_HEIGHT - 30)
        button.tag = 'right'
        self.buttons.append(button)

    def on_draw(self):
        """renders the screen"""

        arcade.start_render()

        # Draw sprites
        self.scaffold_list.draw()
        self.r_sprite_list.draw()
        self.draw_values()

        # Draw the menu bar
        # arcade.draw_rectangle_filled(INVENTORY_WIDTH / 2,
        #                              SCREEN_HEIGHT,
        #                              INVENTORY_WIDTH,
        #                              self.vh,
        #                              color=arcade.color.OXFORD_BLUE)

        arcade.draw_rectangle_filled(SCREEN_WIDTH / 2,
                                     SCREEN_HEIGHT,
                                     SCREEN_WIDTH,
                                     self.vh,
                                     color=arcade.color.OXFORD_BLUE)

        arcade.draw_text('Molecule Builder', int(SCREEN_WIDTH*0.59), SCREEN_HEIGHT - 50, color=arcade.color.WHITE,
                         font_size=30)

        arcade.draw_text(f'Displaying: R{ord(self.tag[0].lower()) - 96}', 10, SCREEN_HEIGHT - self.vh * 0.5 - 20,
                         color=arcade.color.OXFORD_BLUE, font_size=11)

        instructions = ['Welcome to the Drug Discovery Game. Above you can see the starting scaffold',
                        'with the vectors marked by starred numbers. Select r groups from the scrol-',
                        'lable inventory on the left by double clicking to add to the scaffold. You can ',
                        'filter the r groups (descending) by clicking the filter buttons at the top. To see',
                        'the different sets of r groups available for each vector, click the arrows. ',
                        'Change views by using the right and left keys on the keyboard.']
        inst = "Welcome to the Drug Discovery Game. Above you can see the starting scaffold with the vectors marked by starred numbers. Select r groups from the scrollable inventory on the left by double clicking to add to the scaffold. You can filter the r groups (descending) by clicking the filter buttons at the top. To see the different sets of r groups available for each vector, click the arrows. Change views by using the right and left keys on the keyboard."
        instructions = textwrap.fill(inst, 78)
        instructions = instructions.split(sep='\n')
        for i, t in enumerate(instructions):
            arcade.draw_text(t, INVENTORY_WIDTH + 15, SCREEN_HEIGHT / 5 - (i + 1) * 20, color=arcade.color.OXFORD_BLUE)

        # Delineate boundaries
        arcade.draw_line(INVENTORY_WIDTH, SCREEN_HEIGHT, INVENTORY_WIDTH, 0, arcade.color.OXFORD_BLUE, 5)
        arcade.draw_line(INVENTORY_WIDTH, SCREEN_HEIGHT / 5, SCREEN_WIDTH, SCREEN_HEIGHT / 5,
                         arcade.color.OXFORD_BLUE, 5)

        # Draw the filters
        self.filter_sprite_list.draw()
        for f in self.filter_sprite_list:
            if f.tag == self.feature:
                f._set_color(arcade.color.RED)

        # Draw the navigation buttons
        self.buttons.draw()

        # Draw hover text
        if self.display_hover:
            self.draw_hover()

    def on_mouse_press(self, x, y, button, key_modifiers):
        """ Called when the user presses a mouse button. """

        # Find what the user has clicked on
        clicked = arcade.get_sprites_at_point((x, y), self.filter_sprite_list)
        if len(clicked) > 0:
            # Change the clicked filter button
            [f._set_color(arcade.color.WHITE) for f in self.filter_sprite_list]
            feature = clicked[-1]
            feature._set_color(arcade.color.RED)  # turn filter red
            self.feature = feature.tag

            # Sort the r_sprites
            self.desc_df.sort_values(feature.tag, inplace=True, ascending=False)

            # redraw reordered sprites
            self.setup_sprites(self.tag, str(feature.tag))
            self.r_sprite_list.draw()

        # Pick the r sprite and shade it
        for r in self.r_sprite_list:
            r._set_alpha(255)  # removes the shade from any other r_sprites

        if arcade.get_sprites_at_point((x, y), self.r_sprite_list):
            # If already selected, then add to scaffold
            if arcade.get_sprites_at_point((x, y), self.r_sprite_list)[-1] == self.picked_r:
                self.update_lead()
                self.setup()
                self.on_draw()
            self.picked_r = arcade.get_sprites_at_point((x, y), self.r_sprite_list)[-1]  # pick the top sprite
            self.picked_r.smiles = self.desc_df.loc[
                self.desc_df[self.tag] == self.picked_r.tag, 'mol'].item()  # give smile
            self.picked_r._set_alpha(50)  # shade

        # Change inventory
        clicked = arcade.get_sprites_at_point((x, y), self.buttons)
        if len(clicked) > 0:
            [f._set_color(arcade.color.WHITE) for f in self.buttons]
            if clicked[-1] == self.buttons[1]:  # if right arrow
                if ord(self.tag[0]) - 96 < self.num_vecs:  # not out of range
                    clicked[-1]._set_color(arcade.color.RED)
                    self.setup_sprites(tag=f'{chr(ord(self.tag[0]) + 1)}tag', feat=self.feature)
                else:
                    pass
            elif clicked[-1] == self.buttons[0]:  # if left arrow
                if self.tag[0] == 'a':
                    pass
                else:
                    clicked[-1]._set_color(arcade.color.RED)
                    self.setup_sprites(tag=f'{chr(ord(self.tag[0]) - 1)}tag', feat=self.feature)
            else:
                pass

    def on_mouse_scroll(self, x: int, y: int, scroll_x: int, scroll_y: int):
        """Redraw the sprites lower instead of scrollling"""

        for i, r in enumerate(self.r_sprite_list):
            r.position = (r.position[0], r.position[1] + int(scroll_y * 3))
        self.view_top += scroll_y

    def on_update(self, delta_time: float):
        """Checks to see if the user is hovering over a sprite looking for help"""
        # Specify which sprites have help text
        all_sprite_list = self.combine_sprite_lists([self.filter_sprite_list, self.buttons])

        hovered = arcade.get_sprites_at_point(self.location, all_sprite_list)
        self.display_hover = False
        if len(hovered) == 1:
            if self.hovered != hovered[-1]:  # if hovering over something new
                self.hovered = hovered[-1]  # store the sprite that's being hovered over
                self.hover_time = 0
            else:
                self.hover_time += delta_time
            if self.hover_time > 1:
                self.display_hover = True  # feeds back into on_draw()

    def on_mouse_motion(self, x: float, y: float, dx: float, dy: float):
        """Update mouse location"""
        self.location = (x, y)

    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """
        if symbol == arcade.key.R:
            # Restart
            self.setup()

        if symbol == arcade.key.RIGHT:
            if 0 in self.current_rs:
                print("You haven't selected enough r groups")
            else:
                # navigate to feedback buttons view
                pause = FeedbackView(self)  # passes the current view to FeedbackView for later
                self.window.show_view(pause)  # show the Feedback view


# Run the game loop
def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = MolView()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()


if __name__ == "__main__":
    main()
