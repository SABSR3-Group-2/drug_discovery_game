import arcade
import os
import pandas as pd
import re
from descriptors import get_descriptors
import pyglet

"""
Inventory
"""

# Constants
SCREEN_WIDTH = 500
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Inventory"

# Constants for sprite scaling
MOL_SCALING = 0.5
FILTER_SCALING = 0.3
CURSOR_SCALING = 1

# Calculate all descriptors and store in Dataframe
data = pd.read_csv('data/r_group_decomp.csv')  # read data
cols = [c for c in data.columns if re.match('R*\d', c)]  # get the r group cols
desc = [get_descriptors(x) for x in data[cols[0]].unique()]  # calculate descriptors for all moieties of a given r group
desc_df = pd.DataFrame(desc)  # make dataframe
desc_df.insert(0, 'atag', data['atag'].unique())  # insert tags

_window = pyglet.window.Window


class Inventory(arcade.Window):
    """
    Inventory window class
    """

    def __init__(self):
        # call the parent class and set up the window
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)

        # Make relative units for responsive design
        self.vw = int(SCREEN_WIDTH / 6)  # relative width
        self.vh = int(SCREEN_HEIGHT / 6)  # relative height

        # list to hold the sprites
        self.r_sprite_list = None  # r group sprites
        self.filter_sprite_list = None  # feature filter button sprites

        # Track which feature the r groups are being sorted by
        self.feature = None
        self.filters = ['MW', 'logP', 'TPSA', 'HA', 'h_acc', 'h_don', 'rings']

        # Used to keep track of our scrolling
        self.view_top = SCREEN_HEIGHT
        # self.scrolled = 0
        self.top_bound = 0  # the maximum y value
        self.bottom_bound = 0  # the minimum y value

        # Set background colour
        arcade.set_background_color(arcade.color.WHITE)

        # Sprites for parsing
        self.picked_r_list = []

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
        self.r_sprite_list = arcade.SpriteList(use_spatial_hash=True)

        # Sort r_sprites by specified feature
        desc_df.sort_values(feat, inplace=True, ascending=False)  # sort r groups by specified feature
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
        Set up the game here. Call this function to restart the game:

            Process to setting up the r_sprites
                1. Create SpriteList of Sprites
                2. Create list of coordinates
                3. Assigning coordinates to each sprite

        """
        # Set up r group sprites
        # Create a list of the sprites
        self.r_sprite_list = arcade.SpriteList(use_spatial_hash=True)

        # Read in the sprite .pngs and create sprites
        tag = 'atag'
        self.setup_sprites(tag=tag, feat='MW')

        # Set up feature filter buttons
        # Read in filter sprites
        self.filter_sprite_list = arcade.SpriteList(use_spatial_hash=False)
        for i, f in enumerate(self.filters):  # create and position a button for each filter
            filter_sprite = arcade.Sprite(f'Images/filter_pngs/{f}.png', FILTER_SCALING)
            filter_sprite.position = (self.vw * (0.7 * i + 1), self.view_top - 30)
            filter_sprite.filter = f
            self.filter_sprite_list.append(filter_sprite)

    def on_draw(self):
        """Render the screen"""

        # Clear the screen to the background colour
        arcade.start_render()

        # Draw r_sprites
        self.r_sprite_list.draw()

        # Draw the menu bar
        arcade.draw_rectangle_filled(SCREEN_WIDTH / 2,
                                     SCREEN_HEIGHT,
                                     SCREEN_WIDTH,
                                     self.vh,
                                     color=arcade.color.OXFORD_BLUE)

        arcade.draw_text('press c to confirm selection', 0, SCREEN_HEIGHT - int(self.vh / 2) - 20,
                         color=arcade.color.OXFORD_BLUE)

        # Draw the filters
        self.filter_sprite_list.draw()

    def on_mouse_press(self, x: float, y: float, button: int, modifiers: int):
        """Called when user presses a mouse button. Specifically for sorting the r groups by the specified feature"""

        # Updating the feature button
        # Find what the user has clicked on (y coordinate is calculated dynamically)
        clicked = arcade.get_sprites_at_point((x, y), self.filter_sprite_list)

        # Change the clicked button
        if len(clicked) > 0:
            [f._set_color(arcade.color.WHITE) for f in self.filter_sprite_list]
            feature = clicked[0]
            feature._set_color(arcade.color.DARK_CANDY_APPLE_RED)  # turn filter red

            # Sort the r_sprites
            desc_df.sort_values(feature.filter, inplace=True, ascending=False)

            # redraw reordered sprites
            self.setup_sprites('atag', str(feature.filter))
            self.r_sprite_list.draw()
        else:
            pass

        # Pick 3 r groups
        picked_r = arcade.get_sprites_at_point((x, y), self.r_sprite_list)[0]
        if picked_r not in self.picked_r_list:
            if len(self.picked_r_list) < 3:
                picked_r._set_alpha(50)
                self.picked_r_list.append(picked_r)
            else:
                print('you already have 3 groups selected. Unselect one to choose another')
        else:
            picked_r._set_alpha(255)  # remove colour
            self.picked_r_list.remove(picked_r)

    def on_mouse_scroll(self, x: int, y: int, scroll_x: int, scroll_y: int):
        """Redraw the sprites lower instead of scrollling"""
        for i, r in enumerate(self.r_sprite_list):
            r.position = (r.position[0], r.position[1] + scroll_y)
        self.view_top += scroll_y

    def on_key_press(self, symbol: int, modifiers: int):
        """Close the window and save the selected r groups to a tmp file"""
        if symbol == arcade.key.C:
            tags = [r.tag for r in self.picked_r_list]  # get the tags of the selected sprites
            with open('game_scripts/tmp_r_groups.txt', 'w+') as out_file:
                for t in tags:
                    out_file.write(f'{t}\n')
            global _window
            _window.close(self)

    def return_r_groups(self):
        """Return the picked r group list for use in main script"""
        return self.picked_r_list


def inventory_main():
    """Main inventory method"""
    inventory = Inventory()
    inventory.setup()
    arcade.run()


if __name__ == "__main__":
    inventory_main()
