#!/usr/bin/env python3
# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
import re

file = open("custom_taxonomy_table_mqc.txt", "r")
line = file.read()
file.close()
file = open("custom_taxonomy_table_mqc.txt", "w")
line = line.replace("Scaffold", "ViralOTU")
file.write(line)
file.close()
