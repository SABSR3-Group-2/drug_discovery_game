"""Run a game loop"""
import arcade
from molecule_builder import MolView
from preprocessing import preprocess
import re

# Screen title and size
SCREEN_WIDTH = 1000
SCREEN_HEIGHT = 650
INVENTORY_WIDTH = int(SCREEN_WIDTH / 2.2)
SCREEN_TITLE = "Go Forth and Discover Drugs"

# Constants for sprite scaling
CHARACTER_SCALING = 1
TILE_SCALING = 0.5
MOL_SCALING = 0.5
FILTER_SCALING = 0.3
CURSOR_SCALING = 1


def main():
    """Starts the game"""
    usr_input = input("\n----------------------------------------------------------------------------------------------------------------"
                      "\n Welcome to The Drug Discovery Game"
                      "\n----------------------------------------------------------------------------------------------------------------"
                      "\n[1] To quick start using the default MMP-12 inhibitor data, press ENTER."
                      "\n"
                      "\n[2] If you wish to use alternative data which has already been preprocessed, type the file "
                      "name and press ENTER. \n \te.g. my_data.csv (the file must be in the 'data' directory)."
                      "\n"
                      "\n[3] If your data has not yet been preprocessed, please close the game and "
                      "follow the steps in the README."
                      "\n----------------------------------------------------------------------------------------------------------------\n")
    if len(usr_input) < 3:  # start the game
        print("Game starting with default data...")
        pass

    elif re.findall(r'\.csv', usr_input):  # user has provided raw data, now prompt for the scaffold
        scaffold = input("\nPlease provide SMILE string of your scaffold (e.g. O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)["
                         "*:1]) and then press ENTER\n")

        if len(scaffold) < 3:  # user did not provide scaffold correctly
            raise ValueError('Scaffold entered incorrectly')
        else:
            preprocess(raw_data=usr_input, scaffold=scaffold)
    elif usr_input == 'test':
        preprocess('pickett_data.csv', 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
    else:  # user did not provide raw data correctly
        raise ValueError('The file you provided could not be accessed. Please check the file is in the data directory\n'
                         'and you have entered the name correctly')

    # Start the GUI
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
    start_view = MolView()
    window.show_view(start_view)
    start_view.setup()
    arcade.run()


if __name__ == '__main__':
    main()
