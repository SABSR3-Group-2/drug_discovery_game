import arcade

# Screen title and size
SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650

class EndView(arcade.View):
    """Class for generating the view at the end of the game"""

    def __init__(self, mol_view=None):

        # call the parent class and set up the window
        super().__init__()
        self.mol_view = mol_view

        arcade.set_background_color(arcade.color.OXFORD_BLUE)

    def setup(self):
        """This function sets up the screen"""
        pass

    def on_draw(self):
        """Draw the screen"""
        arcade.start_render()

        # Write game over message
        arcade.draw_text("Game Over", SCREEN_WIDTH/2, SCREEN_HEIGHT/2, font_size=70, color=arcade.color.WHITE_SMOKE)