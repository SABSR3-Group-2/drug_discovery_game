"""
Analysis view
"""

import arcade
from end_game_screen import EndView

SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Analysis"

class AnalysisView(arcade.View):
    """
    Analysis class
    """

    def __init__(self, feedback_view=None):
        super().__init__()
        self.feedback_view = feedback_view
        self.end_button_list = None

        arcade.set_background_color(arcade.color.WHITE)

        self.setup()

    def setup(self):
        """
        This function sets up the view, call it to restart.
        """
        # create end button
        self.end_button_list = arcade.SpriteList()
        end_button = arcade.Sprite(f'Images/button_pngs/end_game_blue.png', 0.5)
        end_button.position = SCREEN_WIDTH - 50, 30
        end_button.name = 'end'
        self.end_button_list.append(end_button)

    def on_draw(self):
        """Renders the screen"""
        arcade.start_render()

        # draw the end button
        self.end_button_list.draw()

    def on_mouse_press(self, x, y, button, modifiers):
        """
        Called when the user presses a mouse button. Used for determining what happens
        when the user clicks on a button.
        """

        # identifies what button the user clicks on
        clicked = arcade.get_sprites_at_point((x, y), self.end_button_list)
        if len(clicked) > 0:  # checks a button has been clicked
            if clicked[0].name == 'end':
                # if end button clicked, csv file created and window closed
                self.feedback_view.final_df.to_csv('data/results.csv', index=False)
                end_view = EndView(self.feedback_view.mol_view)  # create end view and pass mol builder view
                self.window.show_view(end_view)

    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """
        # df containing results printed when user presses SPACE
        if symbol == arcade.key.SPACE:
            print(self.feedback_view.final_df)

        if symbol == arcade.key.LEFT:
            self.window.show_view(self.feedback_view)
            arcade.set_background_color(arcade.color.OXFORD_BLUE)


def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = AnalysisView()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()


if __name__ == "__main__":
    main()
