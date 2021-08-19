import os
import arcade
from combine import MolChoose
from descriptors import get_descriptors, lipinski
from filters import run_filters
from rdkit import Chem
import global_vars
from analysis import AnalysisView

"""
Feedback
"""

# Constants
SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
SCREEN_TITLE = "Feedback"

# button names (and costs/duration for assays)
ASSAYS = {'pic50': {'cost': 70, 'duration': 0.5, 'name': 'pIC50'},
          'cl_mouse': {'cost': 7000, 'duration': 3, 'name': 'Clearance (mouse)'},
          'cl_human': {'cost': 9000, 'duration': 3.5, 'name': 'Clearance (humans)'},
          'logd': {'cost': 1000, 'duration': 1.5, 'name': 'LogD'},
          'pampa': {'cost': 700, 'duration': 1, 'name': 'PAMPA'}}
ACTIONS = ['run_assays', 'clear_choices']
CALCULATIONS = ['calculate_descriptors', 'run_filters']
DESC_NAMES = {'MW': 'Molecular weight',
              'logP': 'LogP',
              'TPSA': 'TPSA',
              'HA': 'Heavy atom count',
              'h_acc': 'Number of HBAs',
              'h_don': 'Number of HBDs',
              'rings': 'Number of rings'}

class Button(arcade.Sprite):
    """Sprite button class"""

    def __init__(self, mol, button, scale=1):
        # hold the button name and image
        self.button = button
        self.image_file_name = os.path.join('Images', 'button_pngs', f'{self.button}.png')

        # get molecule information
        self.chosen_mol = mol

        # call the parent class
        super().__init__(self.image_file_name, scale)

    def get_result(self):
        """
        :Returns the assay result when an assay button is clicked.

        :return: assay result for the chosen molecule
        :rtype: string
        """
        # retrieves the appropriate column name
        if self.button == 'cl_human':
            col = 'clearance_human'
        elif self.button == 'cl_mouse':
            col = 'clearance_mouse'
        else:
            col = self.button
        result = self.chosen_mol.at[0, col]
        return str(result)

    def get_cost(self):
        """
        :Returns the cost of the assay when the assay button is clicked

        :return: assay cost
        :rtype: string
        """
        cost = ASSAYS[self.button]['cost']
        return cost

    def get_duration(self):
        """
        :Returns the duration of the assay when the assay button is clicked

        :return: assay duration
        :rtype: string
        """
        duration = ASSAYS[self.button]['duration']
        return duration

    def get_desc(self):
        """
        :Returns a dictionary of descriptors calculated using the descriptors.py script

        :return: calculated descriptors
        :rtype: string
        """
        descriptors = get_descriptors(self.chosen_mol.at[0, 'mol'])
        descriptors.pop('mol')
        for key, val in descriptors.items():  # round to 1 dp
            if key not in ['HA', 'h_acc', 'h_don', 'rings']:
                descriptors[key] = round(float(val), 1)
        return descriptors

    def run_filt(self):
        """
        : Runs the filters in filter.py (PAINS, NIH, BRENK, ZINC) file

        :return: whether the molecule passes the filter and describes any violations
        :rtype: string
        """
        filter_res = run_filters(Chem.MolFromSmiles(self.chosen_mol.at[0, 'mol']))
        return filter_res


class FeedbackView(arcade.View):
    """
    Main application class
    """
    def __init__(self, mol_view=None):
        # call the parent class and set up the window
        super().__init__()
        self.mol_view = mol_view

        # list to hold the button sprites
        self.button_list = None

        # list to hold the mol sprite
        self.mol_sprite_list = None

        # tracks the assays chosen, the results and results to print
        self.assay_choices = None
        self.assay_results = None
        self.assay_results_print = None
        self.assay_choices_print = None

        # track the total cost and duration of the assays selected
        self.total_cost = None
        self.total_duration = None

        # stores the descriptor and filter results
        self.descriptor_results = None
        self.filter_results = None

        # stores the path to the font file
        self.font = os.path.join('fonts', 'arial.ttf')

        # sets the background color
        arcade.set_background_color(arcade.color.OXFORD_BLUE)

        # store the R group tags (will be updated by the molecule builder)
        self.tags = ['A01', 'B01']  # initialise
        for i, t in enumerate(self.mol_view.current_rs):
            if t.tag == 0:
                self.tags[i] = 0
            else:
                self.tags[i] = t.tag

        # track which sprite we're near for displaying help
        self.hovered = None
        self.hover_time = 0
        self.location = (0, 0)
        self.display_hover = False

        # stores the molecule info
        self.mol = None
        self.final_df = None

        self.setup()

    def make_coordinates(self, sprite_no):
        """
        Function to make the coordinates for the assay button sprites.

        :param sprite_no: button number (i.e. 1-5)
        :type sprite_no: int

        :return: two numbers representing 2D coordinates for each sprite
        :rtype: int, int
        """
        y_slot = SCREEN_HEIGHT / 10
        x_slot_width = SCREEN_WIDTH / 5
        x_slot = (sprite_no * x_slot_width) - (x_slot_width / 2) + x_slot_width
        return x_slot, y_slot

    def setup(self):
        """
        Function to set up the feedback view. Creates the molecule and button sprites
        and sets their positions.
        """
        self.button_list = arcade.SpriteList()
        self.assay_results = []
        self.assay_choices = []
        self.assay_results_print = []
        self.assay_choices_print = []
        self.total_cost = 0
        self.total_duration = []
        self.descriptor_results = {}
        self.filter_results = {}
        self.mol_sprite_list = arcade.SpriteList()

        # retrieve molecule information using the r group tags
        for tag in self.tags:
            if 'A' in tag:
                atag = tag
            elif 'B' in tag:
                btag = tag
        self.mol = MolChoose(atag, btag, DataSource=os.path.join('data', 'r_group_decomp.csv'))
        self.mol = self.mol.reset_index(drop=True)

        # make the molecule sprite using the saved image
        mol_sprite = arcade.Sprite(os.path.join('Images', 'game_loop_images',
                                                f'scaffold{self.mol_view.round_count}.png'))
        mol_sprite.position = (SCREEN_WIDTH - (SCREEN_WIDTH / 6)), (SCREEN_HEIGHT - 150)
        self.mol_sprite_list.append(mol_sprite)

        # make the assay buttons (at bottom of the screen)
        for i, assay in enumerate(ASSAYS.keys()):
            assay_button = Button(self.mol, assay, 1)
            assay_button.position = self.make_coordinates(i)
            self.button_list.append(assay_button)

        # make the other four buttons (at top of the screen)
        for i, action in enumerate(ACTIONS + CALCULATIONS):
            action_button = Button(self.mol, action, 0.6)
            action_button.position = (i + (i + 1)) / 18 * SCREEN_WIDTH, (SCREEN_HEIGHT - 90)
            self.button_list.append(action_button)

    def split_text(self, n_words, text):
        """
        Split the hover text into multiple lines so it fits on the screen.
        """
        words = text.split()
        split = [words[i:i + n_words] for i in range(0, len(words), n_words)]
        split = [" ".join(lst) for lst in split]
        return split

    def draw_hover(self):
        """
        Draw the hover-over help text for the assay buttons
        """
        # specify the help text for each button
        text_dict = {'pic50': 'pIC50 represents the negative log of the IC50 (half-maximal inhibitory concentration), a measure of potency',
                     'cl_mouse': 'This assay measures the metabolic clearance of the drug in mice',
                     'cl_human': 'This assay measures the metabolic clearance of the drug in humans',
                     'logd': 'LogD represents the distribution coefficient used to measure lipophilicity',
                     'pampa': 'This assay measures the permeability of the compounds (parallel artificial membrane permeability assay)',
                     'calculate_descriptors': 'Calculate molecular properties of the molecule',
                     'run_filters': 'Run substructure filters (PAINS, NIH, BRENK, ZINC) on the molecule to identify compounds that may result in false positives',
                     'run_assays': 'Run selected assays. This will deduct time and money from your overall balance',
                     'clear_choices': 'Clear selected assays'}

        full_text = text_dict[self.hovered.button]  # what text to write out
        lines = self.split_text(5, full_text)  # split the text into multiple lines so it fits screen

        # set the position of the help text
        if self.hovered.button != 'pampa':  # pampa button text drawn towards the left of the button due to position on the screen
            line_locs = [(self.hovered.position[0] + 40, self.hovered.position[1] - 20 * i) for i in range(len(lines))]
        else:
            line_locs = [(self.hovered.position[0] - 150, self.hovered.position[1] - 20 * i) for i in range(len(lines))]

        # create the text sprite
        for i, line in enumerate(lines):
            text_sprite = arcade.draw_text(line, line_locs[i][0], line_locs[i][1], color=arcade.color.BLACK, font_size=10, font_name=self.font)

            # draw yellow background for the text
            width = text_sprite.width
            height = text_sprite.height
            arcade.draw_rectangle_filled(line_locs[i][0] + width * 0.5, line_locs[i][1] + height * 0.5,
                                         width + 10, height + 10,
                                         color=arcade.color.YELLOW)

            # draw the text
            text_sprite.draw()

    def on_draw(self):
        """
        Render the screen
        """
        # clear the screen to the background colour
        arcade.start_render()

        # draw the chosen molecule section
        arcade.draw_rectangle_filled((SCREEN_WIDTH - (SCREEN_WIDTH / 6)),
                                     (SCREEN_HEIGHT - (4 / 5 * SCREEN_HEIGHT) / 2),
                                     (SCREEN_WIDTH / 3),
                                     (4 / 5 * SCREEN_HEIGHT),
                                     color=arcade.color.BLACK)

        arcade.draw_rectangle_filled((SCREEN_WIDTH - (SCREEN_WIDTH / 6)),
                                     (SCREEN_HEIGHT - (4 / 5 * SCREEN_HEIGHT) / 2),
                                     (SCREEN_WIDTH / 3 - 10),
                                     (4 / 5 * SCREEN_HEIGHT - 10),
                                     color=arcade.color.WHITE)

        arcade.draw_text('Chosen molecule',
                         SCREEN_WIDTH - 260,
                         0.93 * SCREEN_HEIGHT,
                         color=arcade.color.BLACK,
                         font_size=15,
                         font_name=self.font,
                         align='center')

        arcade.draw_text(f"Chosen R groups: {self.mol.at[0, 'atag']}, {self.mol.at[0, 'btag']}",
                         4 / 6 * SCREEN_WIDTH + 20,
                         390,
                         font_size=15,
                         font_name=self.font,
                         color=arcade.color.BLACK)

        self.mol_sprite_list.draw()

        # draw text showing total balances
        arcade.draw_text(f"Total balance: ${global_vars.balance}",
                         4 / 6 * SCREEN_WIDTH + 20,
                         360,
                         font_size=15,
                         font_name=self.font,
                         color=arcade.color.BLACK)

        if global_vars.balance <= 0:  # if balance is negative, text appears in red
            arcade.draw_text(f"Total balance: ${global_vars.balance}",
                             4 / 6 * SCREEN_WIDTH + 20,
                             360,
                             font_size=15,
                             font_name=self.font,
                             color=arcade.color.DARK_CANDY_APPLE_RED)

        arcade.draw_text(f"Time remaining: {global_vars.time} weeks",
                         4 / 6 * SCREEN_WIDTH + 20,
                         330,
                         font_size=15,
                         font_name=self.font,
                         color=arcade.color.BLACK)

        if global_vars.time <= 0: # if time remaining is negative, text appears in red
            arcade.draw_text(f"Time remaining: {global_vars.time} weeks",
                            4 / 6 * SCREEN_WIDTH + 20,
                            330,
                            font_size=15,
                            font_name=self.font,
                            color=arcade.color.DARK_CANDY_APPLE_RED)

        # draw the molecule report section
        arcade.draw_rectangle_filled((1 / 3 * SCREEN_WIDTH),
                                     (1 / 2 * SCREEN_HEIGHT),
                                     (2 / 3 * SCREEN_WIDTH),
                                     (3 / 5 * SCREEN_HEIGHT),
                                     color=arcade.color.BLACK)

        arcade.draw_rectangle_filled((1 / 3 * SCREEN_WIDTH),
                                     (1 / 2 * SCREEN_HEIGHT),
                                     (2 / 3 * SCREEN_WIDTH),
                                     (3 / 5 * SCREEN_HEIGHT - 10),
                                     color=arcade.color.WHITE)

        arcade.draw_text('Molecule report',
                         SCREEN_WIDTH / 3 - 100,
                         SCREEN_HEIGHT * 3 / 4 - 10,
                         font_size=18,
                         font_name=self.font,
                         color=arcade.color.BLACK)

        # draw the instructions text
        instructions = ['Welcome to the feedback screen. Here you',
                        'can run assays on your chosen molecule,',
                        'calculate descriptors and run substructure',
                        'filters. Running assays costs time and money,',
                        'which will be deducted from your total bal-',
                        'ance above. Press the R key to see a summ-',
                        'mary of all molecules made so far, or press',
                        'the L key to return to the molecule builder.']

        for i, t in enumerate(instructions):
            arcade.draw_text(t, 4 / 6 * SCREEN_WIDTH + 20, 320 - (i + 1) * 20, color=arcade.color.BLACK, font_name=self.font)

        # draw the assay results
        arcade.draw_text('Assay results:',
                         30,
                         SCREEN_HEIGHT * 7 / 10 - 25,
                         font_size=15,
                         font_name=self.font,
                         color=arcade.color.BLACK)

        for i, (assa, res) in enumerate(zip(self.assay_choices_print, self.assay_results_print)):
            arcade.draw_text(ASSAYS[assa]['name'],
                             30,
                             SCREEN_HEIGHT - 240 - (i * 20),
                             color=arcade.color.BLACK,
                             font_size=10,
                             font_name=self.font)
            arcade.draw_text(res,
                             160,
                             SCREEN_HEIGHT - 240 - (i * 20),
                             color=arcade.color.BLACK,
                             font_size=10,
                             font_name=self.font)

        # draw the top command buttons
        arcade.draw_text('Commands',
                         10,
                         SCREEN_HEIGHT - 50,
                         font_size=15,
                         font_name=self.font,
                         color=arcade.color.WHITE)

        arcade.draw_text('Free calculations',
                         2 / 9 * SCREEN_WIDTH + 10,
                         SCREEN_HEIGHT - 50,
                         font_size=15,
                         font_name=self.font,
                         color=arcade.color.WHITE)

        self.button_list.draw()

        # draw the total cost and duration of the selected assays
        arcade.draw_text('Total cost to run assays',
                         4 / 9 * SCREEN_WIDTH + 10,
                         SCREEN_HEIGHT - 50,
                         font_size=15,
                         font_name=self.font,
                         color=arcade.color.WHITE)

        cost_text = f"Total cost: ${self.total_cost}"
        arcade.draw_text(cost_text,
                         4 / 9 * SCREEN_WIDTH + 10,
                         SCREEN_HEIGHT - 75,
                         color=arcade.color.WHITE,
                         font_size=10,
                         font_name=self.font)

        if self.total_duration == []:
            duration_text = "Total duration: 0 weeks"
        else:  # assumes assays are run in parallel (records the longest assay in the selection)
            duration_text = f"Total duration: {max(self.total_duration)} weeks"
        arcade.draw_text(duration_text,
                         4 / 9 * SCREEN_WIDTH + 10,
                         SCREEN_HEIGHT - 100,
                         color=arcade.color.WHITE,
                         font_size=10,
                         font_name=self.font)

        # draw descriptor calculations
        arcade.draw_text('Descriptors:',
                         SCREEN_WIDTH * 1 / 3 + 10,
                         SCREEN_HEIGHT * 7 / 10 - 25,
                         color=arcade.color.BLACK,
                         font_size=15,
                         font_name=self.font)

        for i, (desc, val) in enumerate(self.descriptor_results.items()):
            arcade.draw_text(DESC_NAMES[desc],
                             SCREEN_WIDTH * 1 / 3 + 10,
                             SCREEN_HEIGHT - 240 - (i * 20),
                             color=arcade.color.BLACK,
                             font_size=10,
                             font_name=self.font)
            arcade.draw_text(str(val),
                             SCREEN_WIDTH * 1 / 3 + 130,
                             SCREEN_HEIGHT - 240 - (i * 20),
                             color=arcade.color.BLACK,
                             font_size=10,
                             font_name=self.font)

        # draw filter results
        arcade.draw_text('Filters:',
                         30,
                         SCREEN_HEIGHT * 3 / 7 - 30,
                         color=arcade.color.BLACK,
                         font_size=15,
                         font_name=self.font)

        for i, (filt, val) in enumerate(self.filter_results.items()):
            if isinstance(val, list):
                val = ', '.join(val)
                val = val.replace('_', ' ').replace('[', '').replace(']', '')
            arcade.draw_text(filt,
                             30,
                             SCREEN_HEIGHT / 2 - 100 - i * 20,
                             color=arcade.color.BLACK,
                             font_size=10,
                             font_name=self.font)
            arcade.draw_text(str(val),
                             160,
                             SCREEN_HEIGHT / 2 - 100 - i * 20,
                             color=arcade.color.BLACK,
                             font_size=10,
                             font_name=self.font)

        # give result of Lipinski's ro5
        if self.descriptor_results != {}:
            ro5_v, ro5_res = lipinski(self.descriptor_results)
            lipinski_text = f"Molecule {ro5_res} Lipinski filter ({ro5_v} rules broken)"
            arcade.draw_text(lipinski_text,
                             30,
                             SCREEN_HEIGHT / 2 - 100 - 4 * 20,
                             color=arcade.color.BLACK,
                             font_size=10,
                             font_name=self.font)

        # draw hover text
        if self.display_hover:
            self.draw_hover()

    def on_update(self, delta_time: float):
        """
        Checks to see if the user is hovering over a sprite looking for help
        """
        # specify which sprites have help text
        hovered = arcade.get_sprites_at_point(self.location, self.button_list)
        self.display_hover = False
        if len(hovered) == 1:
            if self.hovered != hovered[-1]:  # if hovering over something new
                self.hovered = hovered[-1]  # store the sprite that's being hovered over
                self.hover_time = 0
            else:
                self.hover_time += delta_time
            if self.hover_time > 1:
                self.display_hover = True  # feeds back into on_draw()

    def on_mouse_motion(self, x: float, y: float, dx: float, dy: float):
        """
        Update mouse location
        """
        self.location = (x, y)

    def on_mouse_press(self, x, y, button, modifiers):
        """
        Called when the user presses a mouse button. Used for determining what happens
        when the user clicks on a button.
        """
        # identifies what button the user clicks on
        clicked = arcade.get_sprites_at_point((x, y), self.button_list)
        if len(clicked) > 0:  # checks a button has been clicked
            choice = clicked[0]
            # checks if the button is for an assay
            # the assay name, result, cost and duration are stored

            if choice.button in ASSAYS.keys():
                choice._set_color(arcade.color.YELLOW)  # selected buttons are changed to yellow
                self.assay_choices.append(choice.button)
                self.assay_results.append(choice.get_result())
                self.total_cost += choice.get_cost()
                self.total_duration.append(choice.get_duration())

            # checks if the button is an action button
            elif choice.button in ACTIONS:
                if choice.button == 'run_assays':
                    if self.assay_results == []:
                        pass
                    else:
                        # adds the results to another list to print
                        # changes buttons back to white
                        self.assay_results_print = self.assay_results
                        self.assay_choices_print = self.assay_choices
                        [b._set_color(arcade.color.WHITE) for b in self.button_list]
                        # cost is not deducted if the molecule was not made or assayed
                        if 'Not Made' in self.assay_results_print or 'Not Assayed' in self.assay_results_print:
                            self.total_cost -= ASSAYS['pic50']['cost']
                            self.total_duration.remove(ASSAYS['pic50']['duration'])
                        # costs are subtracted from global variables
                        global_vars.balance -= self.total_cost
                        if len(self.total_duration) >= 1:
                            global_vars.time -= self.total_duration[0]
                        self.assay_results = []
                        self.assay_choices = []
                        self.total_cost = 0
                        self.total_duration = []
                        # append assay data to assay_df
                        # checks if a row for that mol already exists (appends new row if not)
                        if len(self.mol_view.assay_df.loc[
                                (self.mol_view.assay_df['atag'] == self.tags[0]) &
                                (self.mol_view.assay_df['btag'] == self.tags[1])]) == 0:
                            self.mol_view.assay_df = self.mol_view.assay_df.append(
                                {'atag': self.tags[0], 'btag': self.tags[1]}, ignore_index=True)
                        for a, r in zip(self.assay_choices_print, self.assay_results_print):
                            self.mol_view.assay_df.loc[
                                (self.mol_view.assay_df['atag'] == self.tags[0]) &
                                (self.mol_view.assay_df['btag'] == self.tags[1]), a] = r

                elif choice.button == 'clear_choices':
                    # clears the selected assays and recorded data
                    # changes buttons back to white
                    [b._set_color(arcade.color.WHITE) for b in self.button_list]
                    self.assay_results_print = []
                    self.assay_choices_print = []
                    self.assay_results = []
                    self.total_cost = 0
                    self.total_duration = []
                    self.descriptor_results = {}
                    self.filter_results = {}
                    self.assay_choices = []

            # checks if the button is an action button
            elif choice.button in CALCULATIONS:
                if choice.button == 'calculate_descriptors':
                    choice._set_color(arcade.color.YELLOW)
                    self.descriptor_results = choice.get_desc()  # records the descriptor results
                    # append desc data to assay_df
                    # checks if a row for that mol already exists (appends new row if not)
                    if len(self.mol_view.assay_df.loc[
                            (self.mol_view.assay_df['atag'] == self.tags[0]) &
                            (self.mol_view.assay_df['btag'] == self.tags[1])]) == 0:
                        self.mol_view.assay_df = self.mol_view.assay_df.append(
                            {'atag': self.tags[0], 'btag': self.tags[1]}, ignore_index=True)
                    for d, v in zip(self.mol_view.filters, self.descriptor_results.values()):
                        self.mol_view.assay_df.loc[
                            (self.mol_view.assay_df['atag'] == self.tags[0]) &
                            (self.mol_view.assay_df['btag'] == self.tags[1]), d] = v

                elif choice.button == 'run_filters':
                    choice._set_color(arcade.color.YELLOW)
                    self.filter_results = choice.run_filt()  # records the filter results

    def on_key_press(self, key, _modifiers):
        """
        Allow the user to navigate between views using L and R keys
        """
        if key == arcade.key.LEFT:
            if global_vars.balance > 0:
                # navigate back to molecule builder view
                self.window.show_view(self.mol_view)
                arcade.set_background_color(arcade.color.WHITE)
            else:
                print('no more attempts left')

        if key == arcade.key.RIGHT:
            # navigate to view containing analysis (name can be changed)
            self.final_df = self.mol_view.assay_df  # create df that can be passed to AnalysisView
            pause = AnalysisView(self)  # passes the current view to Analysis for later
            self.window.show_view(pause)


def main():
    """ Main method """
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = FeedbackView()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()


if __name__ == "__main__":
    main()
