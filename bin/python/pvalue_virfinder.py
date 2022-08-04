#!/usr/bin/env python3
# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
import os, sys
import pandas as pnd
import re

seqID = sys.argv[1]
inputFile = open('./' + seqID + '_results.txt', 'r') 
exportFile = open('./' + seqID + '_results_tab.txt', 'w')
for line in inputFile:
   new_line = re.sub('^\s+', '', line)
   new_line = re.sub('name ', 'name\t', new_line)
   new_line = re.sub('length\s+', 'length\t', new_line)
   new_line = re.sub('score\s+', 'score\t', new_line)
   new_line = re.sub('^[0-9]*\s+', '', new_line)
   new_line = re.sub('\s{2,}', '\t', new_line)
   new_line = re.sub('(?<=\d) (?=\d)', '\t', new_line)
   exportFile.write(new_line)
inputFile.close()
exportFile.close()
df = pnd.read_csv(
    './' + seqID + '_results_tab.txt', sep='\t', header=0,
    names=['header', 'length', 'score', 'pvalue']
)
f = open('./' + seqID + '_filtered_headers.txt', 'w')
for index, row in df.iterrows():
    if row.pvalue <= 0.05:
        f.write(row.header + '\n')
f.close()
