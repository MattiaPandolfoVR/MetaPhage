# Coded by Gioele Lazzari and Mattia Pandolfo (gioele.lazza@studenti.univr.it, mattia.pandolfo@univr.it)
software = "multiqc_model_editor.py"
version = "0.0.1"

import re, sys

assembler = sys.argv[1]

file = open('count_table_mqc.txt', 'r')
line = file.read()
file.close()

file = open('count_table_mqc.txt', 'w')
line = line.replace("predicted viral sequences", "predicted viral consensus sequences")
line = line.replace("# plot_type: 'table'", "# plot_type: 'table'\n# id: 'count_table'")
line = line.replace("_" + assembler + "_vOTUs_consensus", "")
line = line.replace("CovTools count", "Count table")
line = line.replace("ViralSequence", "ViralOTU")
file.write(line)