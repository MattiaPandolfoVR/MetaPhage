import re

file = open("custom_taxonomy_table_mqc.txt", "r")
line = file.read()
file.close()
file = open("custom_taxonomy_table_mqc.txt", "w")
#line = re.sub(r"_length_\d+", "", line)
line = line.replace("Scaffold", "ViralOTU")
file.write(line)
file.close()
