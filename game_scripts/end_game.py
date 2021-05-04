import os
import arcade
from game_scripts.combine import MolChoose
from game_scripts.descriptors import get_descriptors
from game_scripts.filters import compound_check
from rdkit import Chem
import global_vars


"""
Window for ending the game - displays a comparison of chosen molecule vs final target molecule
"""

# Constants
SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "The Final Selection"

class EndGame(arcade.View):
    """
    Main class for endgame view
    """

    def __init__(self, feedback_view = None):
        #call the parent class and set up the window
        super().__init__()
        self.feedback_view = feedback_view

        # stores the path to the font file
        self.font = os.path.join('fonts', 'arial.ttf')

        # sets the background color
        arcade.set_background_color(arcade.color.OXFORD_BLUE)

    def on_draw(self):
        """Render the screen"""

        # clear the screen to the background colour
        arcade.start_render()

        # draw the chosen molecule section
        arcade.draw_rectangle_filled((SCREEN_WIDTH - (SCREEN_WIDTH / 6)),
                                     (SCREEN_HEIGHT - (4 / 5 * SCREEN_HEIGHT) / 2),
                                     (SCREEN_WIDTH / 3),
                                     (4 / 5 * SCREEN_HEIGHT),
                                     color=arcade.color.BLACK)


def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = EndGame()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()


if __name__ == "__main__":
    main()
