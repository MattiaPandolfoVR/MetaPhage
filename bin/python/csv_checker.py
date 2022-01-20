#!/usr/bin/env python3
from pickle import TRUE
import pandas as pd
import sys
import os
 
def checker(df):
    # checking if first column name is "Sample"
    first_col_name = df.columns[0]
    if first_col_name == "Sample":
        first_col = df['Sample'].tolist()
        # checking if first column element name starts with SRR ..
        if all(col.startswith('SRR') for col in first_col):
            os.rename('metadata_to_check.csv', 'metadata.csv')
            return True
    sys.exit("Your metadata is not correctly formatted. Check the metadata paragraph on the wiki at https://github.com/MattiaPandolfoVR/MetaPhage#input_files")

# reading the csv file using read_csv
df = pd.read_csv('metadata_to_check.csv')
checker(df)
