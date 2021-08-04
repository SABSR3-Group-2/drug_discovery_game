from select_mols import *
import sys
import os
import pandas as pd


def preprocess(raw_data, scaffold):
    """Given a .csv, preprocesses the data to the point where the molecule_builder script can be called.

    :param raw_data: file name of the csv containing the original data. Data should be in a similar format as the
    default data contained in `data/raw_data.csv`.
    :type raw_data: str
    :param scaffold: the core of the molecule as a Smile
    :type scaffold: str (smile)

    :return: Returns none but creates a file containing the processed data
    :rtype: File"""
    raw_data = os.path.join('data', raw_data)
    print('Starting preprocessing. If the file is large this can take some time...')
    processed_data_df = read_mols(raw_data, scaffold)
    print('Finished preprocessing')
    print(processed_data_df.shape)
    print(f'Saving processed data to file: {raw_data}_processed')
    processed_data_df.to_csv(f'{raw_data}_processed.csv', encoding='utf-8', index=False)


if len(sys.argv) < 3:  # throw error
    print("Please make sure you have provided both the raw_data path and the scaffold smile string")
else:
    sys.argv[1].rstrip().lstrip()
    sys.argv[2].rstrip().lstrip()
    preprocess(sys.argv[1], sys.argv[2])



