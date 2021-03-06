"""
Analysis view
"""

import arcade
import pandas as pd
from game_scripts.combine import MolChoose
from rdkit import Chem
from rdkit.Chem import Draw
import os

SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Analysis"

CARD_WIDTH = 130
CARD_HEIGHT = 100

class AnalysisView(arcade.View):
    """
    Analysis class
    """

    def __init__(self, feedback_view=None):
        super().__init__()
        self.feedback_view = feedback_view
        self.end_button_list = None
        self.card_list = None
        self.mat_list = None

        arcade.set_background_color(arcade.color.WHITE)

        self.setup()

    def setup(self):
        """
        This function sets up the view, call it to restart.
        """
        # create end button
        self.end_button_list = arcade.SpriteList()
        end_button = arcade.Sprite(f'Images/button_pngs/end_game_blue.png', 0.5)
        end_button.position = SCREEN_WIDTH / 2, SCREEN_HEIGHT - 35
        end_button.name = 'end'
        self.end_button_list.append(end_button)

        # create cards with molecules on them
        self.card_list = arcade.SpriteList()
        for index, row in self.feedback_view.final_df.iterrows():
            mol_info = MolChoose(row['atag'], row['btag'], DataSource=os.path.join('data', 'r_group_decomp.csv')).reset_index(drop=True)
            mol = Chem.MolFromSmiles(mol_info.at[0, 'mol'])
            Chem.Draw.MolToFile(mol, 'Images/game_loop_images/selected_mol{}.png'.format(index),
                            size=(CARD_WIDTH, CARD_HEIGHT), imageType=None)
            card_sprite = arcade.Sprite('Images/game_loop_images/selected_mol{}.png'.format(index), scale=1)
            card_sprite.position = self.make_coordinates(index + 1)
            self.card_list.append(card_sprite)

        # create blank mats for molecules to sit on
        self.mat_list = arcade.SpriteList()
        for i in range(1, 13):
            mat_sprite = arcade.SpriteSolidColor(width = CARD_WIDTH + 5, height = CARD_HEIGHT + 5, color=arcade.color.LIGHT_BLUE)
            mat_sprite.position = self.make_coordinates(i)
            self.mat_list.append(mat_sprite)

    def make_coordinates(self, index):
        """
        Function to plot the coordinates for the card and mat sprites.
        """
        y_slot_height = (SCREEN_HEIGHT / 6)
        x_slot_1 = CARD_WIDTH / 2 + 10
        x_slot_2 = SCREEN_WIDTH - (CARD_WIDTH / 2 + 10)
        if index <= 6:
            return x_slot_1, y_slot_height * index - y_slot_height / 2
        else:
            return x_slot_2, y_slot_height * (index - 6) - y_slot_height / 2

        return x_slot, y_slot

    def on_draw(self):
        """
        Renders the screen
        """
        arcade.start_render()

        # draw the sprites
        self.end_button_list.draw()
        self.mat_list.draw()
        self.card_list.draw()

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
                self.window.close()

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
