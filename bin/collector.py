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

    results = pnd.DataFrame(data = {
        'Category': [],
        'Value': []})

    #for value in alignments:

    
    
    content = """
# id: "Output from my script'
# section_name: 'Custom data file'
# description: 'This output is described in the file header. Any MultiQC installation will understand it without prior configuration.'
# format: 'csv'
# plot_type: 'bargraph'
# pconfig:
#    id: 'custom_bargraph_w_header'
#    ylab: 'Number of things'
<-- REPLACE -->
"""


    string_eater = io.StringIO()
    datas = results.to_csv( path_or_buf = string_eater, header = False, index = False )
    #content = content.replace("<-- REPLACE -->", string_eater.getvalue())
    content = content.replace("<-- REPLACE -->", "Category_1,374\nCategory_2,229\nCategory_3,39\nCategory_4,253")
    file_report = open("custom_report_mqc.csv", "w")
    file_report.write(content)
    file_report.close()


if __name__ == "__main__":

    parameters = parser.parse_args()

    projectDir = os.path.realpath(__file__).replace("bin/collector.py", "")

    manage(projectDir, parameters.alignments)
