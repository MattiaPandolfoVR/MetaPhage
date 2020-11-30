# Coded by Gioele Lazzari (gioele.lazza@studenti.univr.it)
software = "collector.py"
version = "0.0.1"

import os, argparse, io
import pandas as pnd

parser = argparse.ArgumentParser(
    description = '... '
                  '... '
                  '... .',
    formatter_class = argparse.RawTextHelpFormatter)

options = parser.add_argument_group("Options")

options.add_argument('-v', '--version', action='version', version= software + " v" + version)
# what follow are ALL the parameters passed by nextflow
options.add_argument('-a', '--alignments', nargs="+", default=[], required=True)


def manage(projectDir, alignments):

    results = pnd.DataFrame(data = {'Scaffold': []})
    results = results.set_index('Scaffold')

    for record in alignments:
        f = open(record, "r")
        value = int(f.read())
        f.close()
        sample = record.split("___")[1]
        scaffold = record.split("___")[2].replace(".fasta.txt", "")

        if not scaffold in results.index:
            results.append(pnd.Series(name=scaffold))
        if not sample in results.columns:
            results[sample] = None
        results.at[scaffold, sample] = value
    

    #####################
    # custom_table #######
    #####################
    
    content = """# plot_type: 'table'
# section_name: 'Custom table section'
# description: 'A custom text introduction (a few sentences) for this section.'
<-- REPLACE -->
"""

    string_eater = io.StringIO()
    datas = results.to_csv( path_or_buf = string_eater, sep='\t' )
    content = content.replace("<-- REPLACE -->", string_eater.getvalue())
    file_report = open("custom_table_mqc.txt", "w")
    file_report.write(content)
    file_report.close()


    #####################
    # custom_plot #######
    #####################
    
    content = """id: 'Output from my script'
section_name: 'Custom plot section'
description: 'A custom text introduction (a few sentences) for this section.'
plot_type: 'bargraph'
pconfig:
    id: 'custom_bargraph_w_header'
    ylab: 'Number of things'
data:
"""

    transposed = results.transpose()
    for index, sample in transposed.iterrows():
        content = content + "    " + index + ":\n"
        for scaffold in sample.index:
            content = content + "    " + "    " + scaffold + ": " + str(sample[scaffold]) + "\n"
        
    file_report = open("custom_plot_mqc.yaml", "w")
    file_report.write(content)
    file_report.close()


if __name__ == "__main__":

    parameters = parser.parse_args()

    projectDir = os.path.realpath(__file__).replace("bin/collector.py", "")

    manage(projectDir, parameters.alignments)
