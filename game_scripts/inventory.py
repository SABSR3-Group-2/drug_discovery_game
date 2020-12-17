import arcade
import random

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
        super().__init__(SCREEN_WIDTH,SCREEN_HEIGHT,SCREEN_TITLE)

        # Make relative units for responsive design
        self.vw = SCREEN_WIDTH/6  # relative width
        self.vh = SCREEN_HEIGHT/6  # relative height

        # list to hold the r_group sprites
        self.r_sprite_list = None

        # list to hold the

        # Used to keep track of our scrolling
        self.view_top = SCREEN_HEIGHT

        arcade.set_background_color(arcade.color.WHITE)

    def setup(self):
        """ Set up the game here. Call this function to restart the game"""
        self.r_sprite_list = arcade.SpriteList(use_spatial_hash=True)

        # Create the rows of r_sprites
        """
        To-do:
        1.
            Find out how to do scroll +
            
        3.
            Multiple views (so that can mock up what it would be like to move between main game and inventory)
        
        2. 
            no. of coordinates = len(self.r_sprite_list)
            Start TOP LEFT, move down, then across
            Generate coordinates first, then filter r_sprite_list by feature, then draw sprites at coordiates
        """

        coordinate_list = [[self.vw, self.vh*5],
                           [self.vw, self.vh*3],
                           [self.vw, self.vh],
                           [self.vw*3, self.vh*5],
                           [self.vw*3, self.vh*3],
                           [self.vw*3, self.vh],
                           [self.vw * 5, self.vh * 5],
                           [self.vw * 5, self.vh * 3],
                           [self.vw * 5, self.vh]]

        for coordinate in coordinate_list:
            # Place an r group in the inventory
            r_sprite = arcade.Sprite('../../tutorial/images/items/rgroup A12.png', MOL_SCALING)
            r_sprite.position = coordinate
            self.r_sprite_list.append(r_sprite)

        # Position filter bar

        # Position


    def on_draw(self):
        """Render the screen"""

        # Clear the screen to the background colour
        arcade.start_render()

        # Draw r_sprites
        self.r_sprite_list.draw()

        # Draw the menu bar
        arcade.draw_rectangle_filled(SCREEN_WIDTH/2,
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