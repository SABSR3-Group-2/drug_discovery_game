import arcade
import os
import pandas as pd
import re
from descriptors import get_descriptors

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

# Padding for scrolling
LEFT_VIEWPORT_MARGIN = 250
RIGHT_VIEWPORT_MARGIN = 250
BOTTOM_VIEWPORT_MARGIN = 50
TOP_VIEWPORT_MARGIN = 100

# Calculate all descriptors and store in Dataframe
# data = pd.read_csv('../data/r_group_decomp.csv')  # read data
# cols = [c for c in data.columns if re.match('R*\d', c)]  # get the r group cols
# desc = [get_descriptors(x) for x in data[cols[0]].unique()]  # calculate descriptors for all moieties of a given r group
# desc_df = pd.DataFrame(desc)  # make dataframe
# desc_df.insert(0, 'atag', data['atag'].unique())  # insert tags
#
# desc_df.sort_values('rings', inplace=True, ascending=False)


class MyGame(arcade.Window):
    """
    Main application class
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

        arcade.set_background_color(arcade.color.WHITE)

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
            vh = SCREEN_HEIGHT - (self.vh * i)  # the y coordinate for a single row with 3 columns from the top
            single_row = [[self.vw, vh],
                          [self.vw * 3, vh],
                          [self.vw * 5, vh]]
            for c in single_row:  # append each coordinate at a time
                coordinate_list.append(c)

        # make last row
        if last_row == 1:
            coordinate_list.append([self.vw, SCREEN_HEIGHT - self.vh * (full_rows)])
        else:
            coordinate_list.append([self.vw, SCREEN_HEIGHT - self.vh * (full_rows)])
            coordinate_list.append([self.vw * 3, SCREEN_HEIGHT - self.vh * (full_rows)])

        return coordinate_list

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

        # Read in the sprite .pngs
        for i in range(0, 50):
            r_sprite = arcade.Sprite(f'../Images/r_group_pngs/A{i + 1}.png', MOL_SCALING)
            self.r_sprite_list.append(r_sprite)

        # Create coordinates for the r_sprites
        coordinate_list = self.make_coordinates(len(self.r_sprite_list))

        # Assign coordinates to r_sprites
        for i, sprite in enumerate(self.r_sprite_list):
            sprite.position = coordinate_list[i]

        # Set up feature filter buttons
        # Read in filter sprites
        self.filter_sprite_list = arcade.SpriteList(use_spatial_hash=False)
        for i, f in enumerate(self.filters):  # create and position a button for each filter
            filter_sprite = arcade.Sprite(f'../Images/filter_pngs/{f}.png', FILTER_SCALING)
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
                                     self.view_top,
                                     SCREEN_WIDTH,
                                     self.vh,
                                     color=arcade.color.OXFORD_BLUE)

        # Draw the filters
        self.filter_sprite_list.draw()
        # arcade.draw_text(f'Sort by:\n{self.feature}',
        #                  3,
        #                  self.view_top - 40,
        #                  color=arcade.color.WHITE_SMOKE,
        #                  )
        #
        # arcade.draw_text('MW', self.vw, self.view_top - 30, color=arcade.color.WHITE_SMOKE)
        # arcade.draw_circle_filled(self.vw, self.view_top - 30, 5, arcade.color.WHITE_SMOKE)

    def on_mouse_press(self, x: float, y: float, button: int, modifiers: int):
        """Called when user presses a mouse button"""

        # Find what the user has clicked on
        clicked = arcade.get_sprites_at_point((x, y), self.filter_sprite_list)

        # Change the clicked button
        if len(clicked) > 0:
            [f._set_color(arcade.color.WHITE) for f in self.filter_sprite_list]
            clicked[0]._set_color(arcade.color.DARK_CANDY_APPLE_RED)  # turn filter red

            # Sort the r_sprites
            print(clicked[0].filter)

    def on_mouse_scroll(self, x: int, y: int, scroll_x: int, scroll_y: int):
        """Scroll the screen on mouse scroll"""

        if self.view_top + scroll_y * 3 > SCREEN_HEIGHT:
            pass
        elif self.view_top + scroll_y * 3 - SCREEN_HEIGHT < self.r_sprite_list[-1].position[-1] - 100:
            pass
        else:
            self.view_top = int(self.view_top) + scroll_y * 3

        # do the scrolling
        arcade.set_viewport(0, SCREEN_WIDTH,
                            self.view_top - SCREEN_HEIGHT,
                            self.view_top)
        for i, s in enumerate(self.filter_sprite_list):
            s.position = self.vw * (0.7 * i + 1), self.view_top - 30

def main():
    """Main method"""
    window = MyGame()
    window.setup()
    arcade.run()


if __name__ == "__main__":
    main()
