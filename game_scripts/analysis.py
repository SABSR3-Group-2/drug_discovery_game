"""
Analysis view
"""

import arcade
import pandas as pd
from combine import MolChoose
from rdkit import Chem
from rdkit.Chem import Draw
import os
from end_game_screen import EndView


SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
MENU_WIDTH = int(SCREEN_WIDTH / 3)
SCREEN_TITLE = "Analysis"

CARD_WIDTH = 130
CARD_HEIGHT = 70
MAT_WIDTH = 150
MAT_HEIGHT = 200

class Card():
    """
    Card class to store molecule info and text coordinates
    """
    def __init__(self, card_coordinates, atag, btag, pic50=None, cl_mouse=None, cl_human=None, logd=None, pampa=None):
        self.atag = atag
        self.btag = btag
        self.x_tag = card_coordinates[0] - 30
        self.y_tag = card_coordinates[1] + 35
        self.tag_text = f"{self.atag}, {self.btag}"

        self.pic50 = pic50
        self.x_pic50 = card_coordinates[0] - 70
        self.y_pic50 = card_coordinates[1] - 60
        self.pic50_text = f"pIC50: {self.pic50}"

        self.cl_mouse = cl_mouse
        self.x_cl_mouse = card_coordinates[0] - 70
        self.y_cl_mouse = card_coordinates[1] - 80
        self.cl_mouse_text = f"Cl (mouse): {self.cl_mouse}"

        self.cl_human = cl_human
        self.x_cl_human = card_coordinates[0] - 70
        self.y_cl_human = card_coordinates[1] - 100
        self.cl_human_text = f"Cl (human): {self.cl_human}"

        self.logd = logd
        self.x_logd = card_coordinates[0] - 70
        self.y_logd = card_coordinates[1] - 120
        self.logd_text = f"LogD: {self.logd}"

        self.pampa = pampa
        self.x_pampa = card_coordinates[0] - 70
        self.y_pampa = card_coordinates[1] - 140
        self.pampa_text = f"PAMPA: {self.pampa}"

class AnalysisView(arcade.View):
    """
    Main view class
    """

    def __init__(self, feedback_view=None):
        super().__init__()
        self.feedback_view = feedback_view
        self.button_list = None

        # stores the components of the 'cards'
        self.mol_list = None  # stores the mol sprites (the molecule images on the cards)
        self.mat_list = None  # stores the 'mats' (rectangle representing outside of the card)
        self.text_list = []  # stores the text information for each card (generated with the Card class)

        arcade.set_background_color(arcade.color.WHITE)

        # Make relative units for responsive design
        self.vw = int(MENU_WIDTH / 4)  # relative width
        self.vh = int(SCREEN_HEIGHT / 3)  # relative height

        # Used to keep track of our scrolling
        self.view_top = SCREEN_HEIGHT

        # self.scrolled = 0
        self.top_bound = 0  # the maximum y value
        self.bottom_bound = 0  # the minimum y value

        # get font from font file
        self.font = os.path.join('fonts', 'arial.ttf')

        # stores which molecule card has been chosen
        self.mol_choice = None

        self.setup()

    def setup(self):
        """
        This function sets up the view, call it to restart.
        """
        self.button_list = arcade.SpriteList()
        # create a button to allow user to explore a chosen molecule in the molecule_buidler
        builder_button = arcade.Sprite(f'Images/button_pngs/mol_builder.png', 0.5)
        builder_button.position = (MENU_WIDTH * 1/6), SCREEN_HEIGHT - 70
        builder_button.name = 'builder'
        self.button_list.append(builder_button)

        # create a button to allow the user to run more assays on a chosen molecule in feedback_buttons
        assays_button = arcade.Sprite(f'Images/button_pngs/run_assays.png', 0.5)
        assays_button.position = (MENU_WIDTH * 3/6), SCREEN_HEIGHT - 70
        assays_button.name = 'assays'
        self.button_list.append(assays_button)

        # create the final mol button
        final_button = arcade.Sprite(f'Images/button_pngs/final_choice.png', 0.226)
        final_button.position = (MENU_WIDTH * 5/6), SCREEN_HEIGHT - 70
        final_button.name = 'final'
        self.button_list.append(final_button)

        # create the molecule sprites for the cards
        self.mol_list = arcade.SpriteList()
        for index, row in self.feedback_view.mol_view.assay_df.iterrows():
            # get the molecules that have been built/assayed from the final_df
            mol_info = MolChoose(row['atag'], row['btag'], DataSource=os.path.join('data', 'r_group_decomp.csv')).reset_index(drop=True)
            mol = Chem.MolFromSmiles(mol_info.at[0, 'mol'])
            Chem.Draw.MolToFile(mol, 'Images/game_loop_images/selected_mol{}.png'.format(index),
                            size=(CARD_WIDTH, CARD_HEIGHT), imageType=None)  # save image of the molecule to file
            card_sprite = arcade.Sprite('Images/game_loop_images/selected_mol{}.png'.format(index), scale=1)  # create sprites of the molecules
            self.mol_list.append(card_sprite)

        # create blank 'mats' to represent the outline of the cards
        # the mats will be the clickable item for each card to allow the user to select molecules
        self.mat_list = arcade.SpriteList()
        for (index, row), i in zip(self.feedback_view.mol_view.assay_df.iterrows(), range(len(self.mol_list))):
            mat_sprite = arcade.SpriteSolidColor(width = MAT_WIDTH, height = MAT_HEIGHT, color=arcade.color.LIGHT_BLUE)
            # add tag attributes
            mat_sprite.atag = row['atag']
            mat_sprite.btag = row['btag']
            self.mat_list.append(mat_sprite)

        # create coordinates for the mat sprites
        mat_coordinate_list = self.make_coordinates(len(self.mat_list))
        # use the same coordinates but with a higher y value for the mol sprites
        mol_coordinate_list = []
        for coords in mat_coordinate_list:
            mol_coords = []
            mol_coords.append(coords[0])
            mol_coords.append(coords[1] + 50)
            mol_coordinate_list.append(mol_coords)

        # set the coordinates of the molecule sprites
        for i, sprite in enumerate(self.mol_list):
            sprite.position = mol_coordinate_list[i]

        for i, sprite in enumerate(self.mat_list):
            sprite.position = mat_coordinate_list[i]

        # use the Card class to create objects that store the molecule info and coordinates
        for (index, row), coord in zip(self.feedback_view.mol_view.assay_df.iterrows(), mol_coordinate_list):
            cardtext = Card(coord, row['atag'], row['btag'], row['pic50'], row['cl_mouse'], row['cl_human'], row['logd'], row['pampa'])
            self.text_list.append(cardtext)

    def make_coordinates(self, n_sprites):
        """
        Function to plot the coordinates for the card sprites.
        """
        coordinate_list = []

        n_cards = len(self.feedback_view.mol_view.assay_df)

        full_rows = int(n_sprites / 2)
        last_row = n_sprites % 2

        # make full rows
        for i in range(0, full_rows):
            vh = SCREEN_HEIGHT - (self.vh * (i + 1))  # the y coordinate for a single row with 3 columns from the top
            single_row = [[self.vw, vh],
                          [self.vw * 3, vh]]
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

    def on_draw(self):
        """
        Renders the screen
        """
        arcade.start_render()

        arcade.draw_rectangle_filled(SCREEN_WIDTH * 2/3,
                                SCREEN_HEIGHT / 2,
                                SCREEN_WIDTH * 2/3,
                                SCREEN_HEIGHT,
                                color=arcade.color.OXFORD_BLUE)

        # draw the sprites needed for the cards
        self.mat_list.draw()
        self.mol_list.draw()

        # draw the text that goes on the cards
        for c in self.text_list:
            arcade.draw_text(c.tag_text,
                        c.x_tag,
                        c.y_tag,
                        color=arcade.color.BLACK,
                        font_size=10,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.pic50_text,
                        c.x_pic50,
                        c.y_pic50,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.cl_mouse_text,
                        c.x_cl_mouse,
                        c.y_cl_mouse,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.cl_human_text,
                        c.x_cl_human,
                        c.y_cl_human,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.logd_text,
                        c.x_logd,
                        c.y_logd,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')
            arcade.draw_text(c.pampa_text,
                        c.x_pampa,
                        c.y_pampa,
                        color=arcade.color.BLACK,
                        font_size=8,
                        font_name=self.font,
                        align='center')

        # draw the menu bar
        arcade.draw_rectangle_filled(MENU_WIDTH / 2,
                                     SCREEN_HEIGHT,
                                     MENU_WIDTH,
                                     self.vh,
                                     color=arcade.color.OXFORD_BLUE)

        help_text = ["Select a molecule to investigate further or choose",
                     "your favourite molecule and end the game."]
        for i, line in enumerate(help_text):
            arcade.draw_text(line,
                            15,
                            SCREEN_HEIGHT - 25 - i * 15,
                            color=arcade.color.WHITE,
                            font_size=10,
                            font_name=self.font,
                            align='center')

        # draw the button sprites
        self.button_list.draw()

    def on_mouse_press(self, x, y, button, modifiers):
        """
        Called when the user presses a mouse button. Used for determining what happens
        when the user clicks on a button.
        """
        # check if the user has clicked on a card
        clicked = arcade.get_sprites_at_point((x, y), self.mat_list)
        if len(clicked) > 0:  # checks a button has been clicked
            [b._set_color(arcade.color.WHITE) for b in self.mat_list]
            choice = clicked[0]
            choice._set_color(arcade.color.YELLOW)  # selected buttons are changed to yellow
            self.mol_choice = [choice.atag, choice.btag]  # record the tags of the chosen molecule

        # check if the user has clicked on a button
        clicked = arcade.get_sprites_at_point((x, y), self.button_list)
        if len(clicked) > 0:  # checks a button has been clicked
            # if final choice button clicked, csv file created and end page shown
            if clicked[0].name == 'final':
                self.feedback_view.mol_view.assay_df.to_csv('data/results.csv', index=False)
                pause = EndView(self)  # passes the current view to Analysis for later
                self.window.show_view(pause)

            # if the molecule builder button is clicked, the chosen molecule tags are passed to self.feedback_view.mol_view
            elif clicked[0].name == 'builder':
                if self.mol_choice is not None:
                    sprites = []
                    # get the sprites in the mol builder script that match the tags
                    for sprite in self.feedback_view.mol_view.r_sprite_list:
                        if sprite.tag == self.mol_choice[0]:
                            sprites.append(sprite)
                        if sprite.tag == self.mol_choice[1]:
                            sprites.append(sprite)
                    # update the lead molecule in mol builder with the correct r groups
                    for sprite in sprites:
                        self.feedback_view.mol_view.picked_r = sprite
                        self.feedback_view.mol_view.picked_r.smiles = self.feedback_view.mol_view.desc_df.loc[self.feedback_view.mol_view.desc_df[self.feedback_view.mol_view.tag] == self.feedback_view.mol_view.picked_r.tag, 'mol'].item()  # give smile
                        self.feedback_view.mol_view.picked_r._set_alpha(50)  # shade
                        self.feedback_view.mol_view.update_lead()
                        self.feedback_view.mol_view.setup()
                        self.feedback_view.mol_view.on_draw()

                    # also update the molecule shown in the feedback view
                    for i, t in enumerate(self.mol_choice):
                        self.feedback_view.tags[i] = t
                    self.feedback_view.setup()
                    self.feedback_view.on_draw()

                    # show the mol builder view
                    self.window.show_view(self.feedback_view.mol_view)

            # if the run assays button is clicked, the chosen molecule tags are again passed to both views
            elif clicked[0].name == 'assays':
                if self.mol_choice is not None:
                    sprites = []
                    # get the sprites in the mol builder script that match the tags
                    for sprite in self.feedback_view.mol_view.r_sprite_list:
                        if sprite.tag == self.mol_choice[0]:
                            sprites.append(sprite)
                        if sprite.tag == self.mol_choice[1]:
                            sprites.append(sprite)
                    # update the lead molecule in mol builder with the correct r groups
                    for sprite in sprites:
                        self.feedback_view.mol_view.picked_r = sprite  # pick the top sprite
                        self.feedback_view.mol_view.picked_r.smiles = self.feedback_view.mol_view.desc_df.loc[self.feedback_view.mol_view.desc_df[self.feedback_view.mol_view.tag] == self.feedback_view.mol_view.picked_r.tag, 'mol'].item()  # give smile
                        self.feedback_view.mol_view.picked_r._set_alpha(50)  # shade
                        self.feedback_view.mol_view.update_lead()
                        self.feedback_view.mol_view.setup()
                        self.feedback_view.mol_view.on_draw()
                    
                    # also update the molecule shown in the feedback view
                    for i, t in enumerate(self.mol_choice):
                        self.feedback_view.tags[i] = t
                    self.feedback_view.setup()
                    self.feedback_view.on_draw()

                    # show the feedback view
                    self.window.show_view(self.feedback_view)
                    arcade.set_background_color(arcade.color.OXFORD_BLUE)

    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """
        if symbol == arcade.key.LEFT:
            self.window.show_view(self.feedback_view)
            arcade.set_background_color(arcade.color.OXFORD_BLUE)
    
    def on_mouse_scroll(self, x: int, y: int, scroll_x: int, scroll_y: int):
        """Redraw the sprites lower instead of scrollling"""

        # change the y position of the sprites
        for i, r in enumerate(self.mol_list):
            r.position = (r.position[0], r.position[1] + scroll_y)
        for i, r in enumerate(self.mat_list):
            r.position = (r.position[0], r.position[1] + scroll_y)
        
        # change the y position of the card text
        for c in self.text_list:
            c.y_tag += scroll_y
            c.y_pic50 += scroll_y
            c.y_cl_mouse += scroll_y
            c.y_cl_human += scroll_y
            c.y_logd += scroll_y
            c.y_pampa += scroll_y

        self.view_top += scroll_y

def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = AnalysisView()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()

if __name__ == "__main__":
    main()
