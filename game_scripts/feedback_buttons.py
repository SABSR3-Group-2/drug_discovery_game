import arcade
import pandas as pd
from game_scripts.combine import MolChoose
import arcade.gui
from arcade.gui import UIManager
from game_scripts.descriptors import get_descriptors
from game_scripts.filters import compound_check
from rdkit import Chem

SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Feedback"
# Constants for sizing
BUTTON_SCALE = 1

# How big are the cards?
BUTTON_WIDTH = 190 * BUTTON_SCALE
BUTTON_HEIGHT = 140 * BUTTON_SCALE

test_mol = MolChoose('A01', 'B01', DataSource='data/r_group_decomp.csv')
rgroupdata = pd.read_csv('data/r_group_decomp.csv')


ASSAYS = {
    'pic50': {'cost': 70, 'duration': 0.5},
    'cl_mouse': {'cost': 7000, 'duration': 3},
    'cl_human': {'cost': 9000, 'duration': 3.5},
    'logd': {'cost': 1000, 'duration': 1.5},
    'pampa': {'cost': 700, 'duration': 1}
    }

COMMANDS = ['run_assays', 'clear_choices']

class Button(arcade.Sprite):
    """Button sprite"""

    def __init__(self, assay, scale=1):
        self.assay = assay
        self.image_file_name = f"Images/button_pngs/{self.assay}.png"
        super().__init__(self.image_file_name, scale)
    
    def get_result(self):
        if self.assay == 'cl_human':
            col = 'clearance_human'
        elif self.assay == 'cl_mouse':
            col = 'clearance_mouse'
        else:
            col = self.assay
        
        result = test_mol.at[0,col]
        return str(result)
    
    def get_cost(self):
        cost = ASSAYS[self.assay]['cost']
        return cost

    def get_duration(self):
        duration = ASSAYS[self.assay]['duration']
        return duration

class ActionButton(arcade.Sprite):
    """Button sprite"""

    def __init__(self, action, scale=1):
        self.action = action
        self.image_file_name = f"Images/button_pngs/{self.action}.png"
        super().__init__(self.image_file_name, scale)
    
class MyGame(arcade.Window):
    """
    Main view. Really the only view in this example. """
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)

        self.button_list = None
        self.action_button_list = None
        self.assay_results = None
        self.assay_choices = None
        
        arcade.set_background_color(arcade.color.WHITE)
    
    def make_coordinates(self, sprite_no):
        y_slot = 75
        x_slot_width = SCREEN_WIDTH / 5
        x_slot = (sprite_no * x_slot_width) - (x_slot_width/2) + x_slot_width
        return x_slot, y_slot
    
    def setup(self):
        self.button_list = arcade.SpriteList()
        self.action_button_list = arcade.SpriteList()
        self.assay_results = []
        self.assay_choices = []

        for i, assay in enumerate(ASSAYS.keys()):
            button = Button(assay, BUTTON_SCALE)
            button.position = self.make_coordinates(i)
            self.button_list.append(button)
        
        for i, action in enumerate(COMMANDS):
            action_button = ActionButton(action, BUTTON_SCALE)
            action_button.position = 100, (SCREEN_HEIGHT - 200 - i * 100)
            self.action_button_list.append(action_button)

    def on_draw(self):
        """Render the screen"""

        # Clear the screen to the background colour
        arcade.start_render()

        arcade.draw_rectangle_filled(SCREEN_WIDTH / 2,
                                     SCREEN_HEIGHT / 8,
                                     SCREEN_WIDTH,
                                     170,
                                     color=arcade.color.OXFORD_BLUE)
        
        self.button_list.draw()
        self.action_button_list.draw()
        arcade.draw_text('Assay result', SCREEN_WIDTH-250, SCREEN_HEIGHT-50, color=arcade.color.OXFORD_BLUE, font_size=15)
        arcade.draw_text('Balance', SCREEN_WIDTH-250, SCREEN_HEIGHT-150, color=arcade.color.OXFORD_BLUE, font_size=15)
        arcade.draw_text('Time left', SCREEN_WIDTH-250, SCREEN_HEIGHT-250, color=arcade.color.OXFORD_BLUE, font_size=15)

        res = self.assay_results
        #arcade.draw_text(res, SCREEN_WIDTH-250, SCREEN_HEIGHT-100, color=arcade.color.OXFORD_BLUE, font_size=15)
    
    def on_mouse_press(self, x, y, button, modifiers):
        clicked = arcade.get_sprites_at_point((x, y), self.button_list)
        
        if len(clicked) > 0:
            choice = clicked[0]
            if choice.assay in ASSAYS.keys():
                choice._set_color(arcade.color.DARK_CANDY_APPLE_RED)
                self.assay_results.append(choice.get_result())

def main():
    """ Main method """
    window = MyGame()
    window.setup()
    arcade.run()

if __name__ == "__main__":
    main()

