![Run unit tests](https://github.com/SABSR3-Group-2/drug_discovery_game/workflows/Run%20unit%20tests/badge.svg) [![Documentation Status](https://readthedocs.org/projects/drug-discovery-game/badge/?version=latest)](https://drug-discovery-game.readthedocs.io/en/latest/?badge=latest)

# The Drug Discovery Game

The Drug Discovery Game is an interactive game exploring the process of drug discovery based on actual data from drugs in clinical trials or on the market. Players grow a core fragment from pre-defined vectors, guided by biophysical and pharmacological assay data. The aim is to develop a final compound within a specified budget with a potency close to that of the published compound. This serves as both an educative tool for training medicinal chemists in robust decision making, but also as a tool for training machine learning algorithms in the decisions made by pharmaceutical research scientists. It is developed by researchers from the University of Oxford in collaboration with scientists from Roche.
Latest documentation at: https://drug-discovery-game.readthedocs.io/en/latest/ **Documentation not currently working**
## Installation

Currently the project is not pip-installable, but this functionality will be added in due course.

## Creating a Virtual Environment

## Usage

At the moment, this is very much a minimum viable product (MVP), however it is in constant and continual development, and this repository will be updated frequently.

## Preprocessing New Data

The game is limited to scaffold molecules with two vector positions and your data should resemble `data/picket_data.csv`.
As a minimum, the file should have columns [smiles, Atag, Btag] (not case -sensitive) additionally, data for pic50, clearance_mouse, clearance_human,
logd and pampa are needed for full functionality. Since the names for these values can vary between datasets, please 
navigate to `game_scripts/select_mols.py` and at the bottom of `read_mols()` change the column names (e.g. 'pic50_mmp12') to those that are in your raw_data file.

Ensuring you are in the root directory and replacing the two arguments with the name of your data file and the smile string for your scaffold,
run the preprocessing pipeline by typing a command like the following into your terminal:

`python -m game_scripts/preprocessing.py my_raw_data.csv OC(=Ccc210)`

The preprocessing can take quite some time (the default data takes ~20 mins on a 2015 MacBook Pro). Once it is finished it will
create an output csv in `data/` with a name in the format `my_raw_data_processed.csv`. You can now start the game.

## Starting the game

Ensuring you are in the root directory and have the dependencies installed, type the following into your terminal to start the game:

`python -m game_scripts/gameloop.py`

If you want to use the default MMP12 dataset, simply press `ENTER` as instructed to start the game.

Otherwise, if you are using a different but compatible dataset, you will be prompted to enter the name of the file containing the data 
and subsequnently the smile string for your central scaffold.

A GUI should then open and you can begin playing - good luck!

## Gameplay

The game is built around 4 windows:
<ul>
<li>Molecule building</li>
<li>Virtual assaying</li>
<li>Analysis</li>
<li>Game Over</li>
</ul>

### Molecule Builder

In this window your scaffold will be displayed on the right and you will have an inventory system on the left showing the 
fragments that you can add to your scaffold at the first vector position. After selecting this first position, click the right arrow
at the top of the page to view the second inventory of fragments. You can navigate back and forth by clicking the arrows at the top.
When you have a molecule you are happy with press the `right` arrow key to advance to the next window.

## Contributions

We are always very keen to hear your comments, thoughts and suggestions. Do please open an issue with any feedback.

## License
[MIT](https://choosealicense.com/licenses/mit/)
