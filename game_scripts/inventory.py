import arcade
import random
import os
"""
Inventory
"""

# Constants
SCREEN_WIDTH = 500
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Inventory"

# Constants for sprite scaling
MOL_SCALING = 0.5
CURSOR_SCALING = 1

# Padding for scrolling
LEFT_VIEWPORT_MARGIN = 250
RIGHT_VIEWPORT_MARGIN = 250
BOTTOM_VIEWPORT_MARGIN = 50
TOP_VIEWPORT_MARGIN = 100


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

        # list to hold the r_group sprites
        self.r_sprite_list = None

        # list to hold the

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

        """
        To-do:
        1.
            Find out how to do scroll +
            
        3.
            Multiple views (so that can mock up what it would be like to move between main game and inventory)
            - instead just make this as a window which can be opened by a button click in Matthew's window
        
        2. 
            no. of coordinates = len(self.r_sprite_list)
            Start TOP LEFT, move down, then across
            Generate coordinates first, then filter r_sprite_list by feature, then draw sprites at coordinates
        """

        # Create a list of the sprites
        self.r_sprite_list = arcade.SpriteList(use_spatial_hash=False)
        # mock data for now...
        # r_sprite = arcade.Sprite('../../tutorial/images/items/rgroup A12.png', MOL_SCALING)
        #
        # sorted_r_sprites = []
        # for r in range(0,50):
        #     sorted_r_sprites.append(r_sprite)
        # # sorted_r_sprites = [r_sprite for i in range(0, 50)]  # replace this with sorted list
        # for sprite in sorted_r_sprites:
        #     self.r_sprite_list.append(sprite)

        # Read in the sprite pngs
        for i in range(0,50):
            r_sprite = arcade.Sprite(f'../Images/A{i + 1}.png', MOL_SCALING)
            self.r_sprite_list.append(r_sprite)

        # Create coordinates for the r_sprites
        coordinate_list = self.make_coordinates(len(self.r_sprite_list))

        # Assign coordinates to r_sprites
        for i, sprite in enumerate(self.r_sprite_list):
            sprite.position = coordinate_list[i]

        # # Create coordinates for the r_sprites
        # coordinate_list = [[self.vw, self.vh*5],
        #                    [self.vw, self.vh*3],
        #                    [self.vw, self.vh],
        #                    [self.vw*3, self.vh*5],
        #                    [self.vw*3, self.vh*3],
        #                    [self.vw*3, self.vh],
        #                    [self.vw * 5, self.vh * 5],
        #                    [self.vw * 5, self.vh * 3],
        #                    [self.vw * 5, self.vh]]
        #
        # # Assign coordinates to the r_sprites
        # for coordinate in coordinate_list:
        #     # Place an r group in the inventory
        #     r_sprite = arcade.Sprite('../../tutorial/images/items/rgroup A12.png', MOL_SCALING)
        #     r_sprite.position = coordinate
        #     self.r_sprite_list.append(r_sprite)




    def on_draw(self):
        """Render the screen"""

        # Clear the screen to the background colour
        arcade.start_render()

        # Draw r_sprites
        self.r_sprite_list.draw()

        # Draw the menu bar
        arcade.draw_rectangle_filled(SCREEN_WIDTH / 2,
                                     self.view_top,
                                     SCREEN_WIDTH, self.vh,
                                     color=arcade.color.OXFORD_BLUE)

    def on_mouse_scroll(self, x: int, y: int, scroll_x: int, scroll_y: int):
        """Scroll the screen on mouse scroll"""

        self.view_top = int(self.view_top) + scroll_y * 3

        # do the scrolling
        arcade.set_viewport(0, SCREEN_WIDTH,
                            self.view_top - SCREEN_HEIGHT,
                            self.view_top)


def main():
    """Main method"""
    window = MyGame()
    window.setup()
    arcade.run()


if __name__ == "__main__":
    main()
