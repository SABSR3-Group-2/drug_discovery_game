import arcade
import pandas as pd
from game_scripts.combine import MolChoose
import arcade.gui
from arcade.gui import UIManager
from game_scripts.descriptors import get_descriptors
from game_scripts.filters import compound_check
from rdkit import Chem

# To do
# Tidy up code
# Update filter test
# Fix filter result
# Fix Not Made/Not Assayed exceptions

SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Feedback"
# Constants for sizing
BUTTON_SCALE = 1

# How big are the cards?
BUTTON_WIDTH = 190 * BUTTON_SCALE
BUTTON_HEIGHT = 140 * BUTTON_SCALE

test_mol = MolChoose('A01', 'B05', DataSource='data/r_group_decomp.csv')
test_mol = test_mol.reset_index(drop=True)
rgroupdata = pd.read_csv('data/r_group_decomp.csv')

chosen_mol = Chem.MolFromSmiles(test_mol.at[0,'mol'])
Chem.Draw.MolToFile(chosen_mol, 'Images/button_pngs/chosen_mol.png', size=(300, 300), imageType=None)

ASSAYS = {
    'pic50': {'cost': 70, 'duration': 0.5},
    'cl_mouse': {'cost': 7000, 'duration': 3},
    'cl_human': {'cost': 9000, 'duration': 3.5},
    'logd': {'cost': 1000, 'duration': 1.5},
    'pampa': {'cost': 700, 'duration': 1}
    }

ACTIONS = ['run_assays', 'clear_choices']

CALCULATIONS = ['calculate_descriptors', 'run_filters']

class Button(arcade.Sprite):
    """Creates a button sprite"""

    def __init__(self, button, scale=1):
        self.button = button
        self.image_file_name = f"Images/button_pngs/{self.button}.png"
        super().__init__(self.image_file_name, scale)
    
    def get_result(self):
        if self.button == 'cl_human':
            col = 'clearance_human'
        elif self.button == 'cl_mouse':
            col = 'clearance_mouse'
        else:
            col = self.button
        result = test_mol.at[0,col]
        return str(result)
    
    def get_cost(self):
        cost = ASSAYS[self.button]['cost']
        return cost

    def get_duration(self):
        duration = ASSAYS[self.button]['duration']
        return duration
    
    def get_desc(self):
        descriptors = get_descriptors(test_mol.at[0,'mol'])
        descriptors.pop('mol')
        return descriptors
    
    def run_filt(self):
        filter_res = compound_check(Chem.MolFromSmiles(test_mol.at[0,'mol']))
        return filter_res

class MyGame(arcade.Window):
    """
    Main view. Really the only view in this example. """
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)

        self.button_list = None
        self.assay_choices = None
        self.assay_results = None
        self.assay_results_print = None
        self.total_cost = None
        self.total_duration = None
        self.descriptor_results = None
        self.filter_results = None
        self.mol_sprite_list = None
        self.font = None
        
        arcade.set_background_color(arcade.color.OXFORD_BLUE)
    
    def make_coordinates(self, sprite_no):
        y_slot = SCREEN_HEIGHT/10
        x_slot_width = SCREEN_WIDTH / 5
        x_slot = (sprite_no * x_slot_width) - (x_slot_width/2) + x_slot_width
        return x_slot, y_slot
    
    def setup(self):
        self.button_list = arcade.SpriteList()
        self.assay_choices = []
        self.assay_results_print = []
        self.assay_results = []
        self.total_cost = 0
        self.total_duration = []
        self.descriptor_results = {}
        self.filter_results = []
        self.mol_sprite_list = arcade.SpriteList()
        self.font = 'fonts/arial.ttf'

        # make the mol sprite
        mol_sprite = arcade.Sprite('Images/button_pngs/chosen_mol.png')
        mol_sprite.position = (SCREEN_WIDTH-(SCREEN_WIDTH/6)), (SCREEN_HEIGHT - (4/5 * SCREEN_HEIGHT) / 2)
        self.mol_sprite_list.append(mol_sprite)
        
        # make the assay buttons
        for i, assay in enumerate(ASSAYS.keys()):
            assay_button = Button(assay, BUTTON_SCALE)
            assay_button.position = self.make_coordinates(i)
            self.button_list.append(assay_button)
        
        # make the action buttons
        for i, action in enumerate(ACTIONS + CALCULATIONS):
            action_button = Button(action, 0.6)
            action_button.position = (i + (i+1))/12 * SCREEN_WIDTH, (SCREEN_HEIGHT - 90)
            self.button_list.append(action_button)

    def on_draw(self):
        """Render the screen"""

        # Clear the screen to the background colour
        arcade.start_render()

        # Draw the molecule section
        arcade.draw_rectangle_filled((SCREEN_WIDTH-(SCREEN_WIDTH/6)),
                                     (SCREEN_HEIGHT - (4/5 * SCREEN_HEIGHT) / 2),
                                     (SCREEN_WIDTH / 3),
                                     (4/5 * SCREEN_HEIGHT),
                                     color=arcade.color.BLACK)
        
        arcade.draw_rectangle_filled((SCREEN_WIDTH-(SCREEN_WIDTH/6)),
                                     (SCREEN_HEIGHT - (4/5 * SCREEN_HEIGHT) / 2),
                                     (SCREEN_WIDTH / 3 - 10),
                                     (4/5 * SCREEN_HEIGHT - 10),
                                     color=arcade.color.WHITE)

        arcade.draw_text('Chosen molecule',
                        SCREEN_WIDTH-260,
                        0.9*SCREEN_HEIGHT,
                        color=arcade.color.BLACK,
                        font_size=20,
                        font_name=self.font,
                        align='center')
        
        arcade.draw_text(f"Chosen R groups: {test_mol.at[0,'atag']}, {test_mol.at[0,'btag']}",
                        4/6*SCREEN_WIDTH+20,
                        1/5*SCREEN_HEIGHT+20,
                        font_size=15,
                        font_name=self.font,
                        color=arcade.color.BLACK)
        
        self.mol_sprite_list.draw()

        # Draw the report card
        arcade.draw_rectangle_filled((1/3*SCREEN_WIDTH),
                                     (1/2*SCREEN_HEIGHT),
                                     (2/3*SCREEN_WIDTH),
                                     (3/5*SCREEN_HEIGHT),
                                     color=arcade.color.BLACK)

        arcade.draw_rectangle_filled((1/3*SCREEN_WIDTH),
                                     (1/2*SCREEN_HEIGHT),
                                     (2/3*SCREEN_WIDTH),
                                     (3/5*SCREEN_HEIGHT-10),
                                     color=arcade.color.WHITE)
        
        arcade.draw_text('Molecule report',
                        30,
                        SCREEN_HEIGHT*7/10+10,
                        font_size=20,
                        font_name=self.font,
                        color=arcade.color.BLACK)

        # Draw the top command buttons
        arcade.draw_text('Commands',
                        30,
                        SCREEN_HEIGHT-50,
                        font_size=20,
                        font_name=self.font,
                        color=arcade.color.WHITE)
        
        arcade.draw_text('Free calculations',
                        1/3*SCREEN_WIDTH+30,
                        SCREEN_HEIGHT-50,
                        font_size=20,
                        font_name=self.font,
                        color=arcade.color.WHITE)

        self.button_list.draw()

        # Draw the assay results
        arcade.draw_text('Assay results:',
                        30,
                        SCREEN_HEIGHT*7/10-25,
                        font_size=18,
                        font_name=self.font,
                        color=arcade.color.BLACK)

        for i, (assa, res) in enumerate(zip(self.assay_choices, self.assay_results_print)):
            arcade.draw_text(assa,
                            30,
                            SCREEN_HEIGHT-265-(i*40),
                            color=arcade.color.BLACK,
                            font_size=12,
                            font_name=self.font)
            arcade.draw_text(res,
                            130,
                            SCREEN_HEIGHT-265-(i*40),
                            color=arcade.color.BLACK,
                            font_size=12,
                            font_name=self.font)
        
        cost_text = f"Total cost: ${self.total_cost}"
        arcade.draw_text(cost_text,
                        30,
                        1/5*SCREEN_HEIGHT+40,
                        color=arcade.color.BLACK,
                        font_size=15,
                        font_name=self.font)

        if self.total_duration == []:
            duration_text = "Total duration: 0 weeks"
        else:
            duration_text = f"Total duration: {max(self.total_duration)} weeks"
        arcade.draw_text(duration_text,
                        30,
                        1/5*SCREEN_HEIGHT+20,
                        color=arcade.color.BLACK,
                        font_size=15,
                        font_name=self.font)

        arcade.draw_text('Descriptors:',
                        SCREEN_WIDTH*1/3 + 10,
                        SCREEN_HEIGHT*7/10+10,
                        color=arcade.color.BLACK,
                        font_size=18,
                        font_name=self.font)
        
        for i, (desc, val) in enumerate(self.descriptor_results.items()):
            arcade.draw_text(desc,
                            SCREEN_WIDTH*1/3 + 10,
                            SCREEN_HEIGHT-210-(i*20),
                            color=arcade.color.BLACK,
                            font_size=12,
                            font_name=self.font)
            arcade.draw_text(str(val),
                            SCREEN_WIDTH*1/3 + 110,
                            SCREEN_HEIGHT-210-(i*20),
                            color=arcade.color.BLACK,
                            font_size=12,
                            font_name=self.font)

        arcade.draw_text('Filters',
                        SCREEN_WIDTH*1/3 + 10,
                        SCREEN_HEIGHT*3/7,
                        color=arcade.color.BLACK,
                        font_size=18,
                        font_name=self.font)
        
        for i, filt in enumerate(self.filter_results):
            arcade.draw_text(filt,
                            SCREEN_WIDTH*1/3 + 10,
                            SCREEN_HEIGHT/2-70 - i*20,
                            color=arcade.color.BLACK,
                            font_size=12,
                            font_name=self.font)

    def on_mouse_press(self, x, y, button, modifiers):
        clicked = arcade.get_sprites_at_point((x, y), self.button_list)
        if len(clicked) > 0:
            choice = clicked[0]
            if choice.button in ASSAYS.keys():
                choice._set_color(arcade.color.YELLOW)
                self.assay_results.append(choice.get_result())
                self.total_cost += choice.get_cost()
                self.total_duration.append(choice.get_duration())
                self.assay_choices.append(choice.button)

            elif choice.button in ACTIONS:
                if choice.button == 'run_assays':
                    if self.assay_results == []:
                        pass
                    else:
                        self.assay_results_print = self.assay_results
                        [b._set_color(arcade.color.WHITE) for b in self.button_list]
                elif choice.button == 'clear_choices':
                    [b._set_color(arcade.color.WHITE) for b in self.button_list]
                    self.assay_results_print = []
                    self.assay_results = []
                    self.total_cost = 0
                    self.total_duration = []
                    self.descriptor_results = {}
                    self.filter_results = []
                    self.assay_choices = []

            elif choice.button in CALCULATIONS:
                if choice.button == 'calculate_descriptors':
                    choice._set_color(arcade.color.YELLOW)
                    self.descriptor_results = choice.get_desc()
                elif choice.button == 'run_filters':
                    choice._set_color(arcade.color.YELLOW)
                    self.filter_results = choice.run_filt()

def main():
    """ Main method """
    window = MyGame()
    window.setup()
    arcade.run()

if __name__ == "__main__":
    main()

