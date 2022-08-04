#!/usr/bin/env python3
# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
from pickle import TRUE
import pandas as pd
import sys
import os
 
def checker(df):
    # checking if first column name is "Sample"
    first_col_name = df.columns[0]
    if first_col_name == "Sample":
        os.rename('metadata_to_check.csv', 'metadata.csv')
        return True
    sys.exit("Your metadata is not correctly formatted: first column must be Sample. Check the metadata paragraph on the wiki at https://github.com/MattiaPandolfoVR/MetaPhage#input_files")

# reading the csv file using read_csv
df = pd.read_csv('metadata_to_check.csv')
checker(df)
