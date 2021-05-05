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
    usr_input = input("\n-------------------------------------------------"
                      "\nWelcome to The Drug Discovery Game."
                      "\n-------------------------------------------------"
                      "\nIf your data is already processed, press ENTER to start the game."
                      "\nOtherwise, type the path to the raw data csv file in the data directory (e.g. my_raw_data.csv)"
                      "and then press ENTER.")
    if len(usr_input) < 3:  # start the game
        pass

    elif re.findall(r'\.csv', usr_input):  # user has provided raw data, now prompt for the scaffold
        scaffold = input("Please provide SMILE string of your scaffold (e.g. O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]) "
                         "and then press ENTER")
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
