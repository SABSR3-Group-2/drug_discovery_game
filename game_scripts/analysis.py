"""
Analysis view
"""

import arcade
<<<<<<< HEAD
from end_game_screen import EndView
=======
import pandas as pd
from game_scripts.combine import MolChoose
from rdkit import Chem
from rdkit.Chem import Draw
import os
>>>>>>> 5d330437e7dbb3a936a42644300667b62586b474

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
    Analysis class
    """

    def __init__(self, feedback_view=None):
        super().__init__()
        self.feedback_view = feedback_view
        self.button_list = None
        self.card_list = None
        self.mat_list = None
        self.text_list = []

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

        self.mol_choice = None

        self.setup()

    def setup(self):
        """
        This function sets up the view, call it to restart.
        """
        # create end button
        self.button_list = arcade.SpriteList()
        end_button = arcade.Sprite(f'Images/button_pngs/end_game_blue.png', 0.5)
        end_button.position = SCREEN_WIDTH - 70, 50
        end_button.name = 'end'
        self.button_list.append(end_button)

        builder_button = arcade.Sprite(f'Images/button_pngs/mol_builder.png', 0.5)
        builder_button.position = 50, SCREEN_HEIGHT - 70
        builder_button.name = 'builder'
        self.button_list.append(builder_button)

        assays_button = arcade.Sprite(f'Images/button_pngs/run_assays.png', 0.5)
        assays_button.position = 150, SCREEN_HEIGHT - 70
        assays_button.name = 'assays'
        self.button_list.append(assays_button)

        # create cards with molecules on them
        self.card_list = arcade.SpriteList()
        for index, row in self.feedback_view.final_df.iterrows():
            # get the molecules that have been selected from the df
            mol_info = MolChoose(row['atag'], row['btag'], DataSource=os.path.join('data', 'r_group_decomp.csv')).reset_index(drop=True)
            mol = Chem.MolFromSmiles(mol_info.at[0, 'mol'])
            Chem.Draw.MolToFile(mol, 'Images/game_loop_images/selected_mol{}.png'.format(index),
                            size=(CARD_WIDTH, CARD_HEIGHT), imageType=None)  # save image of the molecule to file
            card_sprite = arcade.Sprite('Images/game_loop_images/selected_mol{}.png'.format(index), scale=1)  # create sprites of the molecules
            self.card_list.append(card_sprite)

        # create blank mats for molecules to sit on
        self.mat_list = arcade.SpriteList()
        for (index, row), i in zip(self.feedback_view.final_df.iterrows(), range(len(self.card_list))):
            mat_sprite = arcade.SpriteSolidColor(width = MAT_WIDTH, height = MAT_HEIGHT, color=arcade.color.LIGHT_BLUE)
            mat_sprite.atag = row['atag']
            mat_sprite.btag = row['btag']

            self.mat_list.append(mat_sprite)
        
        # create coordinates for the cards
        coordinate_list = self.make_coordinates(len(self.card_list))
        card_coordinate_list = []
        for coords in coordinate_list:
            new_coords = []
            new_coords.append(coords[0])
            new_coords.append(coords[1] + 50)
            card_coordinate_list.append(new_coords)

        # set the coordinates of the molecule and mat sprites
        for i, sprite in enumerate(self.card_list):
            sprite.position = card_coordinate_list[i]
        
        for i, sprite in enumerate(self.mat_list):
            sprite.position = coordinate_list[i]
        
        # use the Card class to create objects that store the molecule info and coordinates
        for (index, row), coord in zip(self.feedback_view.final_df.iterrows(), card_coordinate_list):
            cardtext = Card(coord, row['atag'], row['btag'], row['pic50'], row['cl_mouse'], row['cl_human'], row['logd'], row['pampa'])
            self.text_list.append(cardtext)
        

    def make_coordinates(self, n_sprites):
        """
        Function to plot the coordinates for the card sprites.
        """
        coordinate_list = []

        n_cards = len(self.feedback_view.final_df)

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

        self.mat_list.draw()
        self.card_list.draw()

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
        
        arcade.draw_text('Select a molecule to investigate further',
                        15,
                        SCREEN_HEIGHT - 25,
                        color=arcade.color.WHITE,
                        font_size=10,
                        font_name=self.font,
                        align='center')

        # draw the sprites
        self.button_list.draw()

    def on_mouse_press(self, x, y, button, modifiers):
        """
        Called when the user presses a mouse button. Used for determining what happens
        when the user clicks on a button.
        """
        # identifies what button the user clicks on
        clicked = arcade.get_sprites_at_point((x, y), self.mat_list)
        if len(clicked) > 0:  # checks a button has been clicked
            [b._set_color(arcade.color.WHITE) for b in self.mat_list]
            #self.mol_choice = None
            choice = clicked[0]
            choice._set_color(arcade.color.YELLOW)  # selected buttons are changed to yellow
            self.mol_choice = [choice.atag, choice.btag]  # record mol chosen using tag attributes

        # identifies what button the user clicks on
        clicked = arcade.get_sprites_at_point((x, y), self.button_list)
        if len(clicked) > 0:  # checks a button has been clicked
            if clicked[0].name == 'end':
                # if end button clicked, csv file created and window closed
                self.feedback_view.final_df.to_csv('data/results.csv', index=False)
<<<<<<< HEAD
                end_view = EndView(self.feedback_view.mol_view)  # create end view and pass mol builder view
                self.window.show_view(end_view)
=======
                self.window.close()
            
            elif clicked[0].name == 'builder':
                if self.mol_choice is not None:
                    sprites = []
                    for sprite in self.feedback_view.mol_view.r_sprite_list:
                        if sprite.tag == self.mol_choice[0]:
                            sprites.append(sprite)
                        if sprite.tag == self.mol_choice[1]:
                            sprites.append(sprite)
                    for sprite in sprites:
                        self.feedback_view.mol_view.picked_r = sprite  # pick the top sprite
                        self.feedback_view.mol_view.picked_r.smiles = self.feedback_view.mol_view.desc_df.loc[self.feedback_view.mol_view.desc_df[self.feedback_view.mol_view.tag] == self.feedback_view.mol_view.picked_r.tag, 'mol'].item()  # give smile
                        self.feedback_view.mol_view.picked_r._set_alpha(50)  # shade
                        self.feedback_view.mol_view.update_lead()
                        self.feedback_view.mol_view.setup()
                        self.feedback_view.mol_view.on_draw()
                    
                    for i, t in enumerate(self.mol_choice):
                        self.feedback_view.tags[i] = t
                    self.feedback_view.setup()
                    self.feedback_view.on_draw()

                    self.window.show_view(self.feedback_view.mol_view)
            
            elif clicked[0].name == 'assays':
                if self.mol_choice is not None:
                    sprites = []
                    for sprite in self.feedback_view.mol_view.r_sprite_list:
                        if sprite.tag == self.mol_choice[0]:
                            sprites.append(sprite)
                        if sprite.tag == self.mol_choice[1]:
                            sprites.append(sprite)
                    for sprite in sprites:
                        self.feedback_view.mol_view.picked_r = sprite  # pick the top sprite
                        self.feedback_view.mol_view.picked_r.smiles = self.feedback_view.mol_view.desc_df.loc[self.feedback_view.mol_view.desc_df[self.feedback_view.mol_view.tag] == self.feedback_view.mol_view.picked_r.tag, 'mol'].item()  # give smile
                        self.feedback_view.mol_view.picked_r._set_alpha(50)  # shade
                        self.feedback_view.mol_view.update_lead()
                        self.feedback_view.mol_view.setup()
                        self.feedback_view.mol_view.on_draw()
                        
                    for i, t in enumerate(self.mol_choice):
                        self.feedback_view.tags[i] = t
                    self.feedback_view.setup()
                    self.feedback_view.on_draw()
                    self.window.show_view(self.feedback_view)
                    arcade.set_background_color(arcade.color.OXFORD_BLUE)
>>>>>>> 5d330437e7dbb3a936a42644300667b62586b474

    def on_key_press(self, symbol: int, modifiers: int):
        """ User presses key """
        # df containing results printed when user presses SPACE
        if symbol == arcade.key.SPACE:  # debugging code
            print(self.feedback_view.final_df)

        if symbol == arcade.key.LEFT:
            self.window.show_view(self.feedback_view)
            arcade.set_background_color(arcade.color.OXFORD_BLUE)
        
        if symbol == arcade.key.RIGHT:  # debugging code
            print(self.mol_choice[0], self.mol_choice[1])
    
    def on_mouse_scroll(self, x: int, y: int, scroll_x: int, scroll_y: int):
        """Redraw the sprites lower instead of scrollling"""

        for i, r in enumerate(self.card_list):
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
