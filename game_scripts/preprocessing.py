from select_mols import *
import os


def preprocess(raw_data, scaffold):
    """Given a .csv, preprocesses the data to the point where the molecule_builder script can be called."""
    raw_data = os.path.join('data', raw_data)
    print('Starting preprocessing. If the file is large this can take some time...')
    processed_data_df = read_mols(raw_data, scaffold)
    print('Finished preprocessing')
    print(processed_data_df.shape)
