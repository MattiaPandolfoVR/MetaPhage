#!/usr/bin/env python3
# Coded by Gioele Lazzari (gioele.lazzari@univr.it)
software = "graphanalyzer.py"
version = "1.2.1"


# import system libraries:
import sys, os, io, argparse, logging, textwrap
from operator import itemgetter
from copy import deepcopy
import pickle
import statistics as stats


# define a function to communicate with the console with colors: 
def consoleout(lvl, msg):
    colors = {
        'HEADER' : '\033[95m',
        'OKBLUE' : '\033[94m',
        'OKCYAN' : '\033[96m',
        'OKGREEN' : '\033[92m',
        'WARNING' : '\033[93m',
        'FAIL' : '\033[91m',
        'ENDC' : '\033[0m',
        'BOLD' : '\033[1m',
        'UNDERLINE' : '\033[4m'
    }
    if   lvl == 'error': 
        msg = colors['FAIL'] + 'ERROR: ' + colors['ENDC'] + msg + '\n'
        sys.stderr.write(msg)
        sys.stderr.flush()
        sys.exit(1)
    elif lvl == 'warning': 
        msg = colors['WARNING'] + 'WARNING: ' + colors['ENDC'] + msg + '\n'
        sys.stderr.write(msg)
        sys.stderr.flush()
    elif lvl == 'okay': 
        msg = colors['OKGREEN'] + 'OKAY: ' + colors['ENDC'] + msg + '\n'
        sys.stderr.write(msg)
        sys.stderr.flush()
    else: 
        consoleout('error', 'consoleout() called with wrong lvl.')
    


# set the argparser:
parser = argparse.ArgumentParser(description = "This script automatically parse vConTACT2 outputs when using INPHARED as database.",
                                 allow_abbrev=False)
parser.add_argument('-v', '--version', 
                    action='version', 
                    version= software + ' v' + version)
parser.add_argument('--graph',  
                    metavar='FILEPATH', 
                    required=True, 
                    type=str,
                    help='c1.ntw output file from vConTACT2')
parser.add_argument('--csv',    
                    metavar='FILEPATH', 
                    required=True,  
                    type=str,
                    help='genome_by_genome_overview.csv output file from vConTACT2')
parser.add_argument('--metas',  
                    metavar='FILEPATH', 
                    required=True,  
                    type=str,
                    help='data_excluding_refseq.tsv file from INPHARED')
parser.add_argument('-o', '--output', 
                    metavar='PATH', 
                    default='./',  
                    type=str,
                    help='path to the output directory')
parser.add_argument('-s', '--suffix', 
                    metavar='STRING', 
                    default='sample',  
                    type=str,
                    help='suffix to append to every file produced in the output directory')
parser.add_argument('--preprocessed', 
                    metavar='STRING', 
                    default=None,  
                    type=str,
                    help='preprocessed genome_by_genome_overview.csv that already includes info from data_excluding_refseq.tsv (to skip the first step of the script)')


# load the external libraries:
try:
    import pandas as pnd 
except ImportError:
    consoleout('error', "The pandas library was not found.")
try:
    import networkx as net # net.spring_layout() could require also scipy.
except ImportError:
    consoleout('error' , "The networkx or scipy library was not found.")
try: 
    from networkx.drawing.layout import spring_layout
    from networkx.drawing.nx_agraph import graphviz_layout
    # This import requires:
    # Mac: conda install graphviz==2.42.3=h055b950_2 -c conda-forge
    # Mac: conda install pygraphviz==1.6=py37hd4be98e_1 -c conda-forge
except ImportError:
    consoleout('error', "The pygraphviz library was not found.") 
try:
    import holoviews as hv
except ImportError:
    consoleout('error', "The holoviews library was not found.")
try:
    import hvplot.networkx as hvnx
    import hvplot.pandas # hvPlot dynamically adds the Pandas .hvplot() method.
    import hvplot # for the save() method.
    from bokeh.resources import INLINE # for viewing html pages also without an internet connection.
except ImportError:
    consoleout('error', "The hvplot or bokeh library was not found.")
try:
    import panel as pnl
    pnl.extension() # before displaying anything with Panel it is always necessary to load the Panel extension.
except ImportError:
    consoleout('error', "The panel library was not found.")





def clusterExtractor(graph, csv_edit, output_path, string_suffix):

    # prepare results dataframe:
    results = pnd.DataFrame(data = {
        'Scaffold': [],
        'Closer': [], # this is the putative species. 
        'Accession': [],
        'Status': [],
        'VC': [],
        'Level': [],
        'Weight': [],
        'Host': [],
        # Below the putative taxonomy table exluding the species:
        'BaltimoreGroup' : [],
        'Realm' : [],
        'Kingdom' : [],
        'Phylum' : [],
        'Class' : [],
        'Order' : [],
        'Family' : [],
        'Subfamily' : [],
        'Genus' : []
        })


    # Function to append the taxonomy result of a vOTU into the passed results dataframe.
    # "closer" is the accession to search for in the 'Genome' column of "csv_edit" used as database.
    # "scaffold" is a string passaed like "vOTU_123".
    # "vc" is the viral subcluster as determined by vConTACT2. 
    def insertReference(csv_edit, results, scaffold, closer, weight, level, status, vc):

        # extract taxonomy
        baltim = "n.a."
        realm = "n.a."
        kingdom = "n.a."
        phylum = "n.a."
        classs = "n.a."
        order = "n.a."
        family = "n.a."
        subfamily = "n.a."
        genus = "n.a."
        species = "n.a."

        # extract other infos
        host = "n.a."

        # "G": not in the graph.
        # "F": present in the graph but not assigned.
        if level != "F" and level != "G":
            
            baltim     = csv_edit.loc[csv_edit['Genome'] == closer, 'BaltimoreGroup'].values[0]
            realm      = csv_edit.loc[csv_edit['Genome'] == closer, 'Realm'].values[0]
            kingdom    = csv_edit.loc[csv_edit['Genome'] == closer, 'Kingdom'].values[0]
            phylum     = csv_edit.loc[csv_edit['Genome'] == closer, 'Phylum'].values[0]
            classs     = csv_edit.loc[csv_edit['Genome'] == closer, 'Class'].values[0]
            order      = csv_edit.loc[csv_edit['Genome'] == closer, 'Order'].values[0]
            family     = csv_edit.loc[csv_edit['Genome'] == closer, 'Family'].values[0]
            subfamily  = csv_edit.loc[csv_edit['Genome'] == closer, 'Subfamily'].values[0]
            genus      = csv_edit.loc[csv_edit['Genome'] == closer, 'Genus'].values[0]
            species    = csv_edit.loc[csv_edit['Genome'] == closer, 'Species'].values[0]
            
            host = csv_edit.loc[csv_edit['Genome'] == closer, 'Host'].values[0]
        

        # format weight (cloud be float or "n.a."):
        if type(weight) != str:
            weight = str(int(round(weight,0))) 

        # append new row:
        results = results.append({
            'Scaffold': scaffold,
            'Closer': species, # this is the putative species. 
            'Accession': closer,
            'Status': status,
            'VC': vc,
            'Level': level,
            'Weight': weight,
            'Host': host,
            # Below the putative taxonomy table exluding the species:
            'BaltimoreGroup' : baltim,
            'Realm' : realm,
            'Kingdom' : kingdom,
            'Phylum' : phylum,
            'Class' : classs,
            'Order' : order,
            'Family' : family,
            'Subfamily' : subfamily,
            'Genus' : genus
            }, ignore_index = True) 

        return results


    ####################
    # PREPARATORY PART #
    ####################

    # create the .log file
    logging.basicConfig(filename = output_path + 'graphanalyzer_' + string_suffix + '.log', 
                        filemode='w', level = logging.INFO, # level sets the threshold
                        format = '%(asctime)s %(levelname)s: %(message)s',
                        datefmt = '%H:%M:%S') 
    logging.info('Processing of vConTACT2 output is started.\n')

    # less problematic when called as variable:
    csv_edit = csv_edit.rename(columns={"VC Subcluster": "VCSubcluster"})
    csv_edit = csv_edit.rename(columns={"VC Status": "VCStatus"})

    # sobstitute NA with '' to avoid future runtime errors:
    csv_edit = csv_edit.fillna('')

    # make a copy of the original csv_edit:
    csv_edit_copy = csv_edit.copy(deep = True)

    # add a column that is a copy of "Genome": it will be updated at every iteration of the algorithm.
    csv_edit_copy['Genome_editable'] = csv_edit['Genome']

    # make a copy of the original graph: this will contain an EDITABLE FLAG ("assignment") for each scaffold node.
    graph_copy = deepcopy(graph)
    for node in graph_copy.nodes: # create an "assignment" attribute for all nodes.
        net.set_node_attributes(graph_copy, {node: node}, "assignment") # store a copy of node's name (Accession).

    ###################
    # PROCESSING PART #
    ###################

    # extract scaffolds from reference genomes:
    scaffolds_total = csv_edit_copy[csv_edit_copy['Genome_editable'].str.contains('vOTU_')]

    # extract scaffolds contained in the graph (as pnd dataframe):
    scaffolds_ingraph = scaffolds_total.copy(deep = True)
    for row in scaffolds_ingraph.itertuples():
        if graph_copy.has_node(row.Genome) == False : 
            
            # remove row from scaffolds_ingraph if scaffold is not in the graph:
            scaffolds_ingraph = scaffolds_ingraph.drop(scaffolds_ingraph[scaffolds_ingraph.Genome == row.Genome].index)
            
            # In this case "status" should be always "Singleton" (that is: never present in the graph).
            status = row.VCStatus # SINGLETON
            # append dropped row in results as level="G" (meaning: NOT IN GRAPH):
            results = insertReference(csv_edit, results, scaffold=row.Genome, closer="n.a.", weight="n.a.", level="G", status=status, vc="n.a.")

    # log first stats:                
    logging.info("Viral scaffolds in total: %i" % len(scaffolds_total))
    logging.info("Viral scaffolds in the graph: %i\n" % len(scaffolds_ingraph))


    # to count how many iterations the algorithm does:
    counter_iterations = 0

    # to count how many new assignments at every iteration:
    counter_new = -1 # -1 is just to start the algorithm.
    while(counter_new != 0): # agorithm stops at the cycle with 0 new assigments.

        counter_new = 0 # reset the counter at every iteration.
        counter_iterations += 1

        logging.info("##############################################")
        logging.info("################ Iteration %i #################" % counter_iterations)
        logging.info("##############################################\n")

        
        # for each viral scaffold in the graph:
        for row in scaffolds_ingraph.itertuples():
            
            # skip already assigned scaffolds:
            # scaffolds_ingraph header:    Genome      Genome_editable
            # assigned scaffold:          {vOTU_...}  {Staphylococcus...}
            # NOT assigned scaffold:      {vOTU_...}  {vOTU_...}
            if not "vOTU_" in row.Genome_editable:
                continue
            
            # extract scaffold name and VCStatus:
            scaffold = row.Genome
            status = row.VCStatus
            logging.info("Scaffold: %s" % (scaffold))
            logging.info("VCStatus: %s" % (status))
              
            # COMPUTE NEIGHBORS (that are: the directly connected nodes):
            neighbors_list = list(graph_copy.neighbors(scaffold))
            logging.info("Neighbors: %i" % len(neighbors_list))
            

            if (status == "Clustered" or status == "Outlier" or status == "Clustered/Singleton" or "Overlap" in status):
                
                vc = "" # scaffold's viral cluster or subcluster
                vc_ori = row.VCSubcluster
                sameclustered_list = [] # other genomes/scaffolds in the same cluster/subcluster
                
                # EXTRACT SAMECLUSTERED LIST:
                logging.info("Extracting sameclustered (sc)...")

                if status == "Outlier": # not clustered but connected

                    vc = "n.a."
                    sameclustered_list = []

                elif status == "Clustered":

                    vc = row.VCSubcluster
                    # extract all scaffolds clustered in the same subcluster (hereafter: the sameclustered)
                    # note for the future: maybe it's better to replace '.str.contains(vc)' with '== vc'
                    # sameclustered = csv_edit_copy[csv_edit_copy['VCSubcluster'].str.contains(vc)]
                    sameclustered = csv_edit_copy[csv_edit_copy['VCSubcluster'] == vc]
                    # remove current scaffold's row
                    sameclustered = sameclustered.drop(sameclustered[sameclustered.Genome == scaffold].index)
                    sameclustered_list = sameclustered["Genome"].tolist()

                elif status == "Clustered/Singleton": # they have a subcluster of their own
                    
                    # extract all the possible subclusteres within the cluster
                    vc = "VC_" + row.VCSubcluster.split("_")[1] + "_"
                    sameclustered = csv_edit_copy[csv_edit_copy['VCSubcluster'].str.contains(vc)]
                    # remove current scaffold's row
                    sameclustered = sameclustered.drop(sameclustered[sameclustered.Genome == scaffold].index)
                    sameclustered_list = sameclustered["Genome"].tolist()

                elif "Overlap" in status: # "Overlap" belong to n cluster

                    # extract the n cluster and consider every of thier subcluster
                    vc_list = status.replace("Overlap (", "").replace(")", "").split("/")
                    sameclustered = pnd.DataFrame() # empty df
                    for vc in vc_list: # extract the rows and glue them the the previous
                        sameclustered = sameclustered.append(csv_edit_copy[csv_edit_copy['VCSubcluster'].str.contains(vc + "_")])
                    # remove current scaffold's row
                    sameclustered = sameclustered.drop(sameclustered[sameclustered.Genome == scaffold].index)
                    sameclustered_list = sameclustered["Genome"].tolist()

                else:
                    consoleout("error", "Strange level found during the extraction of the cluster.")

                # log some stats
                logging.info("Status: %s" % status)
                logging.info("Number of sc: %i" % len(sameclustered_list))
                # check if all sameclustered are contained in the neighbours
                logging.info("Are all sc connected to scaffold?: %s" % set(sameclustered_list).issubset(neighbors_list))

                
                # GET SAMECLUSTERED-CONNECTED LIST ORDERED BY WEIGHT
                # associate every CONNECTED sameclustered with its weight
                logging.info("Extracting sc-connected (scc)...")
                connected_list = [] # this will be a list of DICT:
                # {"scaffold": ___, "weight": ___, "assignment": ___}
                for node in sameclustered_list:
                    try: 
                        w = graph_copy.get_edge_data(scaffold, node)["weight"]
                        assignment = graph_copy.nodes[node]['assignment']
                        connected_list.append({"scaffold": node, "weight": w, "assignment": assignment})
                    except: # not every pair of nodes in the same cluster are directly connected in the graph
                        continue
                # ordinate the list by weight (itemgetter is from the operator library)
                connected_list = sorted(connected_list, key=itemgetter('weight'), reverse=True) # reverse for descending order


                # EXTRACT 'CLOSER_NODE' (= the haviest connected-sameclustered)
                # note that closer_node could be another scaffold!
                try: 
                    closer_node = connected_list[0]["scaffold"] # to store the node connected with the heavier edge
                    closer_node_w = connected_list[0]["weight"] # to store the heavier edge
                except: # not always there is a best_bet: for example "Outlier" are not clustered
                    closer_node = "n.a."
                    closer_node_w = 0.0 
                logging.info("closer_node: %s" % closer_node)
                logging.info("closer_node_w: %f" % closer_node_w)

                
                # 'CLOSER_REF' SEARCH
                closer_ref = "n.a." # store the reference genome connected with the heavier edge
                closer_ref_w = 0.0 # store the heavier egde that connect to a reference genome
                level = "n.a." # keep track of the cardinality (order of the list)

                # iterate the list until the first reference genome is reached
                for i in range(len(connected_list)):
                    # if it's connected a true reference genome or a scaffold assigned in the last iteration:
                    if  connected_list[i]["scaffold"].find("vOTU_") == -1 or connected_list[i]["assignment"].find("vOTU_") == -1:
                        
                        if connected_list[i]["scaffold"].find("vOTU_") == -1:
                            closer_ref = connected_list[i]["scaffold"]
                            closer_ref_w = connected_list[i]["weight"]
                            level = str(counter_iterations) + "C" + str(i + 1) # first level is 1, not 0

                        elif connected_list[i]["assignment"].find("vOTU_") == -1 and counter_iterations == 1:
                            continue # scaffold assigned are not considerated in the first iteration. This will skip the "break"

                        elif connected_list[i]["assignment"].find("vOTU_") == -1:
                            closer_ref = connected_list[i]["assignment"]
                            closer_ref_w = connected_list[i]["weight"]
                            level = str(counter_iterations) + "C" + str(i + 1) # first level is 1, not 0
                        
                        # break the for-loop, because the first reference genome was reached
                        break
                        

                # if a reference genome was found within the connected-sameclustered:            
                if closer_ref != "n.a.": 
                    counter_new += 1

                    logging.info("VCSubcluster: %s" % (vc))
                    logging.info("closer_ref: %s" % closer_ref)
                    logging.info("closer_ref_w: %f" % closer_ref_w) 
                    logging.info("level: %s" % level)          

                    # update results:
                    if status == "Clustered/Singleton":
                        # 'vc_ori' instead of 'vc', beacuse 'vc' changed during the algorithm
                        results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc_ori)
                    elif "Overlap" in status:
                        # 'n.a.' instead of 'vc', beacuse 'vc' changed during the algorithm
                        results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc="n.a.")
                    else:
                        results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc)
                    # upadate starting table and graph
                    scaffolds_ingraph.loc[scaffolds_ingraph["Genome"] == scaffold, ["Genome_editable"]] = closer_ref
                    graph_copy.nodes[scaffold]['assignment'] = closer_ref
                
                
                else: # RETRY ITERATING THE WHOLE NEGHBORS 
                    
                    # associate every neighbors with its weight
                    connected_list = [] # this will be a list of DICT:
                    # {"scaffold": ___, "weight": ___, "assignment": ___}

                    for node in neighbors_list:
                        w = graph_copy.get_edge_data(scaffold, node)["weight"]
                        assignment = graph_copy.nodes[node]['assignment']
                        connected_list.append({"scaffold": node, "weight": w, "assignment": assignment})
                    # ordinate the list by weight (itemgetter is from the operator library)
                    connected_list = sorted(connected_list, key=itemgetter('weight'), reverse=True) # reverse for descending order
                    
                    
                    # 'CLOSER_REF' SEARCH
                    # iterate the list until the first reference genome is reached
                    for i in range(len(connected_list)):
                        if  connected_list[i]["scaffold"].find("vOTU_") == -1 or connected_list[i]["assignment"].find("vOTU_") == -1:

                            if connected_list[i]["scaffold"].find("vOTU_") == -1:
                                closer_ref = connected_list[i]["scaffold"]
                                closer_ref_w = connected_list[i]["weight"]
                                level = str(counter_iterations) + "N" + str(i + 1) # first level is 1, not 0

                            elif connected_list[i]["assignment"].find("vOTU_") == -1 and counter_iterations == 1:
                                continue # scaffold assigned are not considerated in the first iteration. This will skip the "break"
                            
                            elif connected_list[i]["assignment"].find("vOTU_") == -1:
                                closer_ref = connected_list[i]["assignment"]
                                closer_ref_w = connected_list[i]["weight"]
                                level = str(counter_iterations) + "N" + str(i + 1) # first level is 1, not 0
                            
                            # break the for-loop, because the first reference genome was reached
                            break


                    # if a reference genome was found within the whole neighbors:                                 
                    if closer_ref != "n.a.": 
                        counter_new += 1
                        
                        logging.info("VCSubcluster: %s" % (vc))
                        logging.info("closer_ref: %s" % closer_ref)
                        logging.info("closer_ref_w: %f" % closer_ref_w) 
                        logging.info("level: %s" % level)          

                        # update results
                        if status == "Clustered/Singleton":
                            # 'vc_ori' instead of 'vc', beacuse 'vc' changed during the algorithm
                            results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc_ori)
                        elif "Overlap" in status:
                            # 'n.a.' instead of 'vc', beacuse 'vc' changed during the algorithm
                            results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc="n.a.")
                        else:
                            results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc)
                        # upadate starting table and graph
                        scaffolds_ingraph.loc[scaffolds_ingraph["Genome"] == scaffold, ["Genome_editable"]] = closer_ref
                        graph_copy.nodes[scaffold]['assignment'] = closer_ref


                    else: # this means: NO reference in sameclustered AND NO reference in neighbors
                        logging.warning("Iteration %i: NO reference in connected sameclustered AND NO reference in neighbors." % counter_iterations)
               
            else: # "Singleton" never present in the graph
                logging.consoleout('error', "Found a strange VCStatus! Maybe it's a Singleton? Singletons should not be present inside the graph.")   

            logging.info("##############################################\n")


    # add the remaining scaffolds (level="F": present in the graph but still not assigned)
    for row in scaffolds_ingraph.itertuples():
        if "vOTU_" in row.Genome_editable:
            # From what we've seen, status for "F" can be:
            # Clustered, Outlier, Overlap. But NOT "Singleton" (that is: never present in the graph).
            # Please note: also "Clustered/Singleton" is valid VCStatus, though never found in 'F' here till now.
            status = row.VCStatus # should NOT be a "Singleton". 
            vc = row.VCSubcluster
            if status == "Outlier" or "Overlap" in status: 
                vc = "n.a."
            # At this point, row.Genome_editable should still be == to row.Genome.
            results = insertReference(csv_edit, results, scaffold=row.Genome_editable, closer="n.a.", weight="n.a.", level="F", status=status, vc=vc)



    logging.info('Processing of vConTACT2 output is ended.') 
    

    # order results by scaffold name
    results = results.sort_values(by=['Scaffold'])


    # Save the original version in Excel (just for help in the developing):
    results.to_excel(output_path + "results_vcontact2_NOCUT_" + string_suffix + ".xlsx")


    # Cut taxonomy at a level dependent of the "confidence" of the assignment:
    def taxonomyCutter(results):
        # create a true copy of the dataframe
        res = deepcopy(results)
        for index, row in res.iterrows():
            # We decided to stop at:
            # Genus for Subclusters;
            # Subfamily (if exists, otherwise Family) for Clusters;
            # Order for 'N' (Outliers).
            if "C" in row.Level and row.Status == "Clustered":
                continue # the best case: keep everything
            elif "C" in row.Level and (row.Status == "Clustered/Singleton" or "Overlap" in row.Status):
                res.at[index,'Genus'] = "O" # Please note: "O" stands for "omitted".
            elif "N" in row.Level :
                res.at[index,'Genus'] = "O" 
                res.at[index,'Subfamily'] = "O"
                res.at[index,'Family'] = "O"
        return res
    results = taxonomyCutter(results)
    consoleout("okay", "Finished to refine taxonomy levels.")


    # Now we want to save results dataframe into 4 different formats:
    # Formats are: 1) html widget; 2) csv file; 3) Excel file; 4) MultiQC custom table.

    # 1) Save results as an html widget table:
    html_widget_table = results.hvplot.table(width = 800, height = 200)

    # 2) save results as a csv file: 
    results.to_csv(output_path + "results_vcontact2_" + string_suffix + ".csv")

    # 3) save results as a csv file: 
    results.to_excel(output_path + "results_vcontact2_" + string_suffix + ".xlsx")

    # 4) Save results as a MultiQC custom table:
    # First, add the custom MultiQC header to the final string "content":
    content = textwrap.dedent("""
    # id: 'taxotab'
    # plot_type: 'table'
    # section_name: 'vConTACT2 extender table'
    # description: 'Automatic processing of vConTACT2 outputs `c1.ntw` and `genome_by_genome_overview.csv`.'
    <-- REPLACE -->
    """)
    # Crate a StringIO object to flush results.to_csv() into:
    string_eater = io.StringIO()
    datas = results.to_csv(path_or_buf = string_eater, sep='\t', index=False)
    # Get out the results.to_csv() as string, and insert it into "content":
    content = content.replace("<-- REPLACE -->", string_eater.getvalue())
    # Write the whole "content" string on file.
    file_report = open(output_path + "custom_taxonomy_table_" + parameters.suffix + "_mqc.txt", "w")
    file_report.write(content) 
    file_report.close()

    return html_widget_table, results



def plotCreatorGraphvizHoloviews(graph, csv_edit, results, output_path, string_suffix, max_weight):

    # less problematic when called as variable
    csv_edit = csv_edit.rename(columns={"VC Subcluster": "VCSubcluster"})
    csv_edit = csv_edit.rename(columns={"VC Status": "VCStatus"})

    # sobstitute NA with '' to avoid future runtime errors
    csv_edit = csv_edit.fillna('')


    # PART 1.
    # Below we add attributes to each node in the graph. Thanks to this, the interactive
    # visualization will be more informative (each attribute is shown when a node is hovered
    # with the mouse).

    # create empty dict:
    attribs = {}

    # for each row in csv_edit, preparate a col and a dict to rename nodes
    for index, row in csv_edit.iterrows():

        # not all references in genome_by_genome_overview.csv are contained in the graph!
        # for genomes/scaffolds not included in graph: ignore.
        if graph.has_node(row.Genome) == False:
            continue

        # don't want to rename scaffold at this point
        if "vOTU_" in row.Genome: # Discriminate accesions from vOTUs.

            # find the corresponding vOTU in the results table:
            matches = results[results['Scaffold'] == row.Genome]
            if len(matches) != 1:
                # there should be only 1 exact match
                consoleout("error", "We have len(matches) != 1 in plotCreatorGraphvizHoloviews().")
            else:
                # get first (and only) match:
                match = matches.iloc[0]

                # take infos from the results dataframe:
                label= {"A0_Type": "Scaffold",
                        "A1_Species/Closer": match.Closer,
                        "A2_Accession": match.Accession,
                        "A3_Status": match.Status,
                        "A4_VC": match.VC,
                        "A5_Level": match.Level,
                        "A6_Weight": match.Weight,
                        "A7_Genus": match.Genus,
                        "A8_Family": match.Family,
                        "A9_Host": match.Host}
                
                # fill the "dict of dicts":
                attribs[row.Genome] = label


        else: # here we have a reference genome:
            
            # take infos from the csv_edit dataframe:
            label= {"A0_Type": "Reference",
                    "A1_Species/Closer": row.Species,
                    "A2_Accession": row.Accession, # commentable because the index is already the accession.
                    "A3_Status": row.VCStatus,
                    "A4_VC": row.VCSubcluster,
                    "A5_Level": "-",
                    "A6_Weight": "-",
                    "A7_Genus": row.Genus,
                    "A8_Family": row.Family,
                    "A9_Host": row.Host}

            # fill the "dict of dicts":
            attribs[row.Genome] = label

    # add attributes to nodes using the dict just prepared:
    net.set_node_attributes(graph, attribs)



    # PART 2.
    # Here we want to draw the interactive plot (whole graph).

    # Below we want to subset results:
    df_results_ingraph = deepcopy(results)    # not G
    df_results_assigned = deepcopy(results)   # not G AND not F
    df_results_confident = deepcopy(results)  # only 1Cx

    # subset results dataframe: obtain scaffolds "assigned" and scaffolds "in-graph".
    # codes: not assigned ("F") or not-in-graph ("G")
    for row in results.itertuples():
        if row.Level == "F" or row.Level == "G": 
            df_results_assigned = df_results_assigned.drop(df_results_assigned[df_results_assigned.Scaffold == row.Scaffold].index)
            df_results_confident = df_results_confident.drop(df_results_confident[df_results_confident.Scaffold == row.Scaffold].index)
            if row.Level == "G":
                df_results_ingraph = df_results_ingraph.drop(df_results_ingraph[df_results_ingraph.Scaffold == row.Scaffold].index)
        elif not "1C" in row.Level:
            df_results_confident = df_results_confident.drop(df_results_confident[df_results_confident.Scaffold == row.Scaffold].index)


    # extract scaffold that are in-graph / that have received a taxanomy / that are confident:
    scaffolds_ingraph = df_results_ingraph['Scaffold'].tolist()
    scaffolds_assigned = df_results_assigned['Scaffold'].tolist()
    scaffolds_confident = df_results_confident['Scaffold'].tolist()
    

    # Extract all nodes connected (directly or indirectly) with "scaffolds_ingraph" nodes.
    # With this strategy there will be much less nodes to render.
    nodes_connected = []
    for element in scaffolds_ingraph:
        nodes_connected = nodes_connected + list(net.node_connected_component(graph, element))
    nodes_connected = list(set(nodes_connected)) # remove duplicates if present.
    
    graph_connected = graph.subgraph(nodes_connected) # create a sugraph for speeding up subsequent computations.
    

    # exclude scaffolds_assigned from the scoffolds_ingraph
    # So the following list will be equal to all the "F" scaffolds (that is: in graph but not assigned).
    scaffolds_UNassigned = list(set(scaffolds_ingraph)-set(scaffolds_assigned)) # only F
    # exclude scaffolds from the reference genomes
    non_scaffolds = list(set(nodes_connected)-set(scaffolds_ingraph)) # only reference genomes


    # separate "assigneds" from "confidents" to avoid conflicts:
    scaffolds_assigned = list(set(scaffolds_assigned) - set(scaffolds_confident)) # not G and not F and not 1CX


    # So at this point we have the following lists:
    # scaffolds_assigned           : not G AND not F AND not 1CX
    # scaffolds_UNassigned         : only F
    # scaffolds_ingraph            : not G
    # scaffolds_confident          : only 1CX
    # non_scaffolds                : reference genomes (only "connected")


    # LABEL CREATION (ony for "confident" scaffolds):
    labels = df_results_confident['Closer'].tolist()
    labdict = {} # create dictionary for labels:
    for i in range(len(scaffolds_confident)): # This works because it's all orderd.
        labdict[scaffolds_confident[i]] = labels[i]


    # calculate position of nodes. Below some parameters you can pass as strings to graphviz.
    # -K: roughly corresponds to an ideal edge length (in inches). 
    # -len: can be used to override K for adjacent nodes.
    # -repulsiveforce: values larger than 1 tend to reduce the warping effect at the expense of less clustering.
    # -overlap: determines if and how node overlaps should be removed. If "true" , overlaps are retained. 
    pos = graphviz_layout(graph_connected, prog="sfdp",  args='-Goverlap=true')
    # note: sfdp requires graphviz built with gts, but this type of build is rarely available in conda,
    # so sfdp will print "Error: remove_overlap: Graphviz not built with triangulation library".
    # Right builds are graphviz==2.42.3=h055b950_2 and pygraphviz==1.6=py37hd4be98e_1 from conda-forge.


    # PLOTTING as described in https://hvplot.holoviz.org/user_guide/NetworkX.html .
    # colors available at https://docs.bokeh.org/en/latest/docs/reference/colors.html .
    # In holoviews and in hvplot, the * operator superimposes the different layers of the plot.
    # First of all, draw all the EDGES with the same style:
    image = hvnx.draw_networkx_edges(graph_connected, pos=pos, 
                                    edge_width = 0.1, alpha = 0.3, 
                                    width = 800, height = 500)
    # then draw the reference genomes:
    image = image * hvnx.draw_networkx_nodes(graph.subgraph(non_scaffolds), pos=pos, 
                                    node_color="aquamarine", node_size=200, alpha=0.3, linewidths=0.0) 
    # draw scaffold not assigned (== only "F"):
    image = image * hvnx.draw_networkx_nodes(graph.subgraph(scaffolds_UNassigned), pos=pos, 
                                    node_color="blue", node_size=200, alpha=0.5, linewidths=0.0) 
    # draw assigned scaffolds (== not G AND not F AND not 1CX):
    image = image * hvnx.draw_networkx_nodes(graph.subgraph(scaffolds_assigned), pos=pos, 
                                    node_color="orange", node_size=200, alpha=0.5, linewidths=0.0)
    # draw "confident" assigned scaffolds with labels (== only 1CX):
    image = image * hvnx.draw_networkx_nodes(graph.subgraph(scaffolds_confident), pos=pos, 
                                    node_color="red", node_size=200, alpha=0.5,  linewidths=0.0, 
                                    labels=labdict, font_size = "8pt")
    

    # save the interactive graph to html before returning it:
    hvnx.save(image, output_path + "graph_layout_" + parameters.suffix + '.html')
    # below we want to edit previous html to fit the MultiQC custom html format.
    # The new edited html will be named "graph_layout_mqc.html".
    file = open(output_path + "graph_layout_" + parameters.suffix + '.html', "r")
    line = file.read()
    file.close() # always close file streams!
    # MultiQC wants some directive as html comment placed above:
    line = textwrap.dedent("""
    <!--
    id: 'graphpane'
    section_name: 'vConTACT2 graph explorer'
    description: 'Interactive graph. Use buttons to interact. Hover nodes to see references.'
    -->
    """) + line # prepend the MultiQC header.
    file = open(output_path +  "graph_layout_" + parameters.suffix + '_mqc.html', "w")
    file.write(line)
    file.close() # always close file streams!

    consoleout("okay", "Finished to generate the whole interactive graph.")




    # PART 3.
    # Here we want to draw a neighbours-based plot for each vOTUs in graph.

    # First of all create the subfolder for storing all the subplots.
    desired_path = output_path + 'single-views_' + string_suffix + '/'
    if (os.path.isdir(desired_path) == False): # if it's not been created yet:
        try:
            os.mkdir(desired_path) # make it!
        except:
            consoleout("error", "Can't create the output sub-folder '%s'. " % desired_path)
    
    # Extract all vOTUs inside the graph: already done in 'scaffolds_ingraph'.

    # For each 'scaffolds_ingraph', compute the nieghtbors:
    for scaffold in scaffolds_ingraph:

        neigh = list(graph.neighbors(scaffold))
        # recreate a subgraph with: nieghbors + scaffold
        # we assume that nieghbors() doesn't already include 'scaffold'.
        neigh.append(scaffold)
        sview_graph = graph.subgraph(neigh)
        

        # the following is just debugging purposes:
        debug_mode = False
        if debug_mode == True:
            if(scaffold == "vOTU_1"): # - Clustered 1C1 ("all to all")
                with open(output_path + 'debugging_vOTU_1.bin', 'wb') as vOTU_1_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_1_file)
            if(scaffold == "vOTU_5"): # - Clustered/Singleton 1C1
                with open(output_path + 'debugging_vOTU_5.bin', 'wb') as vOTU_5_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_5_file)
            if(scaffold == "vOTU_6"): # - Clustered 1C2
                with open(output_path + 'debugging_vOTU_6.bin', 'wb') as vOTU_6_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_6_file)
            if(scaffold == "vOTU_10"): # - Outlier 1N1
                with open(output_path + 'debugging_vOTU_10.bin', 'wb') as vOTU_10_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_10_file)
            if(scaffold == "vOTU_13"): # - Outlier 1N1
                with open(output_path + 'debugging_vOTU_13.bin', 'wb') as vOTU_13_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_13_file)
            if(scaffold == "vOTU_15"): # - Clustered/Singleton 1C2
                with open(output_path + 'debugging_vOTU_15.bin', 'wb') as vOTU_15_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_15_file)
            if(scaffold == "vOTU_16"): # - Outlier 1N1  
                with open(output_path + 'debugging_vOTU_16.bin', 'wb') as vOTU_16_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_16_file)
            if(scaffold == "vOTU_20"): # - Outlier 1N2
                with open(output_path + 'debugging_vOTU_20.bin', 'wb') as vOTU_20_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_20_file)
            if(scaffold == "vOTU_30"): # - Outlier 1N4
                with open(output_path + 'debugging_vOTU_30.bin', 'wb') as vOTU_30_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_30_file)
            if(scaffold == "vOTU_35"): # - Clustered 1C2 and rich subgraph
                with open(output_path + 'debugging_vOTU_35.bin', 'wb') as vOTU_35_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_35_file)
            if(scaffold == "vOTU_42"): # - Clustered (some not-sub-clustered in the same figure)
                with open(output_path + 'debugging_vOTU_42.bin', 'wb') as vOTU_42_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_42_file)
            if(scaffold == "vOTU_51"): # - Outlier 1N1
                with open(output_path + 'debugging_vOTU_51.bin', 'wb') as vOTU_51_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_51_file)
            if(scaffold == "vOTU_63"): # - Outlier 1N2
                with open(output_path + 'debugging_vOTU_63.bin', 'wb') as vOTU_63_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_63_file)
            if(scaffold == "vOTU_64"): # - Clustered 1N1
                with open(output_path + 'debugging_vOTU_64.bin', 'wb') as vOTU_64_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_64_file)
            if(scaffold == "vOTU_71"): # - Overlap 1C1
                with open(output_path + 'debugging_vOTU_71.bin', 'wb') as vOTU_71_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_71_file)
            if(scaffold == "vOTU_74"): # - Overlap 1N1
                with open(output_path + 'debugging_vOTU_74.bin', 'wb') as vOTU_74_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_74_file)
            if(scaffold == "vOTU_82"): # - Outlier 1N3
                with open(output_path + 'debugging_vOTU_82.bin', 'wb') as vOTU_82_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_82_file)
            if(scaffold == "vOTU_99"): # - Clustered 1N3
                with open(output_path + 'debugging_vOTU_99.bin', 'wb') as vOTU_99_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_99_file)
            if(scaffold == "vOTU_100"): # - Overlap 1N3
                with open(output_path + 'debugging_vOTU_100.bin', 'wb') as vOTU_100_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_100_file)
            if(scaffold == "vOTU_109"): # - problems while coloring
                with open(output_path + 'debugging_vOTU_109.bin', 'wb') as vOTU_109_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_109_file)
            if(scaffold == "vOTU_114"): # - Clustered/Singleton 1N4
                with open(output_path + 'debugging_vOTU_114.bin', 'wb') as vOTU_114_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_114_file)
            if(scaffold == "vOTU_118"): # - Clustered 1N2
                with open(output_path + 'debugging_vOTU_118.bin', 'wb') as vOTU_118_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_118_file)
            if(scaffold == "vOTU_122"): # - Clustered 1N2 
                with open(output_path + 'debugging_vOTU_122.bin', 'wb') as vOTU_122_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_122_file)
            if(scaffold == "vOTU_140"): # - Overlap 2N1 (no reference)
                with open(output_path + 'debugging_vOTU_140.bin', 'wb') as vOTU_140_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_140_file)
            if(scaffold == "vOTU_146"): # - Clustered 1N2
                with open(output_path + 'debugging_vOTU_146.bin', 'wb') as vOTU_146_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_146_file)
            if(scaffold == "vOTU_156"): # - Outlier 1N3
                with open(output_path + 'debugging_vOTU_156.bin', 'wb') as vOTU_156_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_156_file)
            if(scaffold == "vOTU_159"): # - Overlap 1N2
                with open(output_path + 'debugging_vOTU_159.bin', 'wb') as vOTU_159_file: # wb = write binary.
                    pickle.dump(sview_graph, vOTU_159_file)

        print("Now generating the " + scaffold + " subgraph.")   
        """
        START paste
        """
        # Extract attributes of the current scaffold (vOTU):
        scaff_Type      = sview_graph.nodes[scaffold]["A0_Type"]
        scaff_Accession = sview_graph.nodes[scaffold]["A2_Accession"]
        scaff_Status    = sview_graph.nodes[scaffold]["A3_Status"]
        scaff_VC        = sview_graph.nodes[scaffold]["A4_VC"]
        scaff_Level     = sview_graph.nodes[scaffold]["A5_Level"]

        if scaff_Status == "Singleton" or scaff_Level == 'G':
            consoleout("error", "There shouldn't be 'G' scaffolds at this point. Contact the developer.")


        # Calculate position of all nodes and edges with spring_layout(): 
        # (APPROXIMATELY, the higher is 'weight' attribute, the shorter is edge length) 
        sview_pos = spring_layout(sview_graph, weight='weight')


        # try to set the theme:
        hv.renderer('bokeh').theme = 'dark_minimal'
        # draw empty graph 1000x650:
        sview_image = hvnx.draw_networkx(sview_graph, pos=sview_pos, 
                            edgelist=[], nodelist=[],  
                            width=1000, height=650)


        # update the attribute 'A6_Weight' for every node:
        for node in sview_graph: # here node is a string !
            attrs = {} # a dict.
            if node == scaffold: attrs = {node: {"A6_Weight": "origin"}}
            else: attrs = {node: {"A6_Weight": str(round(sview_graph[scaffold][node]["weight"],1))}}
            net.set_node_attributes(sview_graph, attrs)


        # Now we want to draw all the edges:
        # Drawing edges one-by-one is too slow and generates too heavy .html files:
        """
        for edge in sview_graph.edges.data("weight"): # (nodeA, nodeB, weight)
            sview_image = sview_image * hvnx.draw_networkx_edges(
                    sview_graph, pos=sview_pos, edgelist=[(edge[0], edge[1])],
                    edge_color='grey', # 'grey' for for testing
                    edge_width= edge[2] / max_weight *3)
        """
        # Using 'cmap' and 'dim' is faster and lighter.
        # For example 90 KB vs 1.3 MB for vOTU_1
        """
        sview_image = sview_image * hvnx.draw_networkx_edges(
                sview_graph, pos=sview_pos,
                edge_color='weight', # See https://hvplot.holoviz.org/user_guide/NetworkX.html 
                # edge_cmap='coolwarm', # See https://matplotlib.org/stable/tutorials/colors/colormaps.html
                edge_cmap='brg',
                # Normalize width to max_weight, and *3 to make thicker lines:
                edge_width=hv.dim('weight') / max_weight *3 )
        """
        # Another strategy could be to pass to 'edge_color' and 'edge_width' of hvnx.draw_networkx_edges
        # an ordered list (color_list, width_list). But we begin creating a list of dict, 
        # because it's sortable based on some key (for future uses). Starting from the list of dict we'll
        # create the necessary lists (color_list, width_list, ...).
        edge_list = [] # list of dict
        for edge in sview_graph.edges.data("weight"): # (nodeA, nodeB, weight)
            curr_dict = {'pairs': None, 'weight': None, 'width': None, 'color': None, 'alpha': None} 
            curr_dict['pairs'] = (edge[0], edge[1]) # 'Tuples' are written with round brackets.
            curr_dict['weight'] = edge[2]
            curr_dict['width'] = edge[2] / max_weight *6 # '*6' as a scaling factor.
            curr_dict['alpha'] = 0.3 if edge[2] / max_weight < 0.3 else 1.0
            # This is a "3-point" gradient, and thus below we use (255*2).
            # The 3 points are: (0,255,255) ; (255,255,0) ; (255,0,0)
            factor = int(round(edge[2] /max_weight * (255*2)))
            if (factor <= 255): curr_dict['color'] = "#%02x%02x%02x" % (factor, 255, 255-factor)
            else: curr_dict['color'] = "#%02x%02x%02x" % (255, 255-(factor-255), 0)
            edge_list.append(curr_dict) # finally append the dict
        # convert the list of dict to the necessary lists (color_list, width_list, ...).
        pairs_list, width_list ,color_list, alpha_list = [], [], [], []
        for edge in edge_list: 
            pairs_list.append(edge['pairs'])
            width_list.append(edge['width'])
            color_list.append(edge['color'])
            alpha_list.append(edge['alpha'])
        # Draw all the edges:
        sview_image = sview_image * hvnx.draw_networkx_edges(
                sview_graph.edge_subgraph(pairs_list), pos=sview_pos, edgelist=pairs_list,
                edge_color=color_list, edge_width=width_list, alpha=alpha_list)



        for node in sview_graph: # here node is a string !

            # Extract attributes of the current node:
            curr_Type        = sview_graph.nodes[node]["A0_Type"]
            curr_Accession   = sview_graph.nodes[node]["A2_Accession"]
            curr_Status      = sview_graph.nodes[node]["A3_Status"]
            curr_VC          = sview_graph.nodes[node]["A4_VC"]
            curr_Level       = sview_graph.nodes[node]["A5_Level"]

            # Properties of this node:
            curr_shape = ""
            curr_size  = 0
            curr_color = "orchid"

            # Distinguish between "Scaffold" and "Reference":
            if curr_Type == "Scaffold":
                curr_shape = "circle"
                curr_size = 400
            elif curr_Type == "Reference": 
                curr_shape = "triangle"
                curr_size = 400
            else:
                consoleout("error", "Strange A0_Type when determining the node's shape.")

            """
            conditions = [  "1C" in scaff_Level,
                            "2C" in scaff_Level and scaff_Status == "Clustered",
                            "1N" in scaff_Level and scaff_Status == "Clustered",
                            "1N" in scaff_Level and "Overlap" in scaff_Status,
                            "2N" in scaff_Level and "Overlap" in scaff_Status,
                            "2N" in scaff_Level and scaff_Status == "Clustered",
                            scaff_Level == "F" and "Overlap" in scaff_Status,
                            scaff_Level == "F" and scaff_Status == "Clustered"
                        ]
            """
            conditions = [True]


            # SCAFFOLD view point.
            valid_scaff_VCs = []
            if any(conditions):
                # Clustered 
                if scaff_Status == "Clustered":
                    valid_scaff_VCs = [scaff_VC]
                # Clustered/Singleton
                elif scaff_Status == "Clustered/Singleton":
                    valid_scaff_VCs = ["VC_" + scaff_VC.split("_")[1] + "_"] 
                # Outlier
                elif scaff_Status == "Outlier": # Outliers are never put inside a VC:
                    valid_scaff_VCs = []
                # Overlap
                elif "Overlap" in scaff_Status:
                    valid_scaff_VCs = scaff_Status.replace("Overlap (", "").replace(")", "").split("/") 

            
            # CURRENT NODE view point:
            valid_curr_VCs = []
            if any(conditions):
                # Clustered, Clustered/Singleton 
                if curr_Status == "Clustered" or curr_Status == "Clustered/Singleton": # VC_z_k
                    if "Overlap" in scaff_Status: # VC_z1, VC_z2, VC_z3
                        valid_curr_VCs = ["VC_" + curr_VC.split("_")[1]] 
                    elif scaff_Status == "Clustered/Singleton": # VC_z_*
                        valid_curr_VCs = ["VC_" + curr_VC.split("_")[1] + "_"]
                    else: 
                        valid_curr_VCs = [curr_VC] # VC_z_k
                # Outlier
                elif curr_Status == "Outlier": # Outliers are never put inside a VC:
                    valid_curr_VCs = [] # Remember: ">>> [] in []" return 'False'
                # Overlap
                elif "Overlap" in curr_Status: # VC_z1, VC_z2, VC_z3
                    valid_curr_VCs = curr_Status.replace("Overlap (", "").replace(")", "").split("/")
                    if scaff_Status == "Clustered" or scaff_Status == "Clustered/Singleton":
                        valid_scaff_VCs = ["VC_" + vc.split("_")[1] for vc in valid_scaff_VCs]
                        

            # color for nodes in the same VC, distinguishing "real" VCs (Clustered) from 
            # "artifical" or "extended" VCs (Overlap and CLustered/Singleton):
            if any(conditions):
                if any(vc in valid_scaff_VCs for vc in valid_curr_VCs):
                    if "Overlap" in scaff_Status or scaff_Status == "Clustered/Singleton":
                        curr_color = "yellow" 
                    else: 
                        curr_color = "darkorange"


            # Add Cluster's nodes when in the Subcluster mode:
            # So we have to consider the 3 cases: 'Clustered', 'Clustered/Singleton', 'Overlap'. 
            if scaff_Status == "Clustered" and (curr_Status == "Clustered" or curr_Status == "Clustered/Singleton"):
                if (("VC_" + scaff_VC.split('_')[1] + "_") in curr_VC) and (curr_VC != scaff_VC): 
                    curr_color = "yellow" 
            elif scaff_Status == "Clustered" and ("Overlap" in curr_Status):
                # Keep in mind that it's a string like "Overlap (VC_4/VC_412/VC_41)""
                if (("VC_" + scaff_VC.split('_')[1] + "/") in curr_Status): 
                    curr_color = "yellow" 
                elif (("VC_" + scaff_VC.split('_')[1] + ")") in curr_Status): 
                    curr_color = "yellow" 


            # understand if this node determines the taxonomy of the current 'scaffold':
            # this if statement picks up only 'References'. So we'll have just 1 magenta triangle.
            if scaff_Accession == curr_Accession and curr_Type == "Reference":
                curr_color = "limegreen"
            
            # check if this is current vOTU:
            if node == scaffold:
                curr_color = "orangered"

            # draw this node:
            sview_image = sview_image * hvnx.draw_networkx_nodes(
                    sview_graph.subgraph([node]), pos=sview_pos, 
                    node_color=curr_color, node_shape=curr_shape, node_size=curr_size, 
                    alpha=1.0, linewidths=1.0)


        # save this interactive subgraph:
        hvnx.save(sview_image, desired_path + scaffold + '.html')

        # add some html tags to help the user:
        file = open(desired_path + scaffold + '.html', "r")
        wholetext = file.read(); file.close() # always close file streams!
        tags = textwrap.dedent("""
        <body><p>Interactive plot generated with <strong>graphanalyzer.py</strong>. Please wait the loading.</p>
        <p>User guide available at <a href="https://www.github.com/lazzarigioele/graphanalyzer/">github.com/lazzarigioele/graphanalyzer</a>.</p>
        <p>Bugs can be reported to <a href= "mailto:gioele.lazzari@univr.com">gioele.lazzari@univr.com</a>.</p>
        """)
        file = open(desired_path + scaffold + '.html', "w")
        file.write(wholetext.replace("<body>", tags))
        file.close() # always close file streams!
        """
        END paste
        """

    consoleout("okay", "Finished to generate the  neighbors-based plot for each vOTUs.")



    return image # this is the whole graph.


 
def fillWithMetas(csv , metas, output_path, string_suffix):

    # vConTACT2 used with the method described in https://github.com/RyanCook94/inphared.pl
    # produces a c1.ntw with GenBank accessions as node labels, and a genome_by_genome_overview.csv
    # with GenBank accessions as genome names. Moreover columns “Genus”, “Family” and “Order” 
    # are filled with "Unassigned" so they need to be filled with the correct taxonomy. 
    # Fortunately, method described at https://github.com/RyanCook94/inphared.pl provides a 
    # metadata table, 26Jan2021_data_excluding_refseq.tsv, containing the correct taxonomy for 
    # every GenBank accession. 

    # The following function update the void columns of genome_by_genome_overview.csv with the
    # taxonomies contained in 26Jan2021_data_excluding_refseq.tsv.

    # make an editable copy:
    csv_edit = csv.copy(deep=True) 

    # rename some columns, so below they could easily be called as Series.name :
    metas = metas.rename(columns={"Sub-family": "Subfamily"})
    metas = metas.rename(columns={"Baltimore Group": "BaltimoreGroup"})

    # Add some new columns:
    csv_edit["Accession"] = csv_edit["Genome"] # copy to allow future edits of "Genome".
    csv_edit["FullClass"] = None # to save the "full classification" as stored in metas.Classification.
    # Host column will be imported from metas, matching the accession number:
    csv_edit["Host"] = None

    # Below we try to reconstruct the full taxonomy. 
    # Genus, Family, and Order are commented as they were already created in the original csv.
    # Genus and Family will anyway be filled with info from metas, matching the accession number. 

    # Taxonomic levels signed with *, are not present in metas tables. We will try to fill
    # these fields decomposing the metas.Classification string, that is a full taxonomy string in ascending order, like eg:
    # "Escherichia phage phiX174 Sinsheimervirus Bullavirinae Microviridae Petitvirales Malgrandaviricetes Phixviricota Sangervirae Monodnaviria Viruses"
    
    csv_edit["BaltimoreGroup"] = None # Not a real taxonomy level - to be filled with info from metas, matching the accession number.
    csv_edit["Realm"] = None # Or "domain" - to be filled with info from metas, matching the accession number.
    csv_edit["Subrealm"] = None # *
    csv_edit["Kingdom"] = None # *
    csv_edit["Subkingdom"] = None # *
    csv_edit["Phylum"] = None # *
    csv_edit["Subphylum"] = None # *
    csv_edit["Class"] = None # *
    csv_edit["Subclass"] = None # *
    # Order # *
    csv_edit["Suborder"] = None # *
    # Family 
    csv_edit["Subfamily"] = None  # to be filled with info from metas, matching the accession number.
    # Genus 
    csv_edit["Species"] = None # to be filled with info from metas, matching the accession number.


    # Viral taxonomy, defined by ICTV, has several taxonomy levels. 
    # Each taxonomic level is charaterized by a particular suffix.
    # The species level has no suffix, can contain more words, and 
    # must not only contain the word virus and the host name. 
    # Below he report all levels-suffix pairs, in form of python dictionary.
    # (level-suffix pairs taken from https://en.wikipedia.org/wiki/Virus_classification)
    lvlsuffix = {"Realm": "viria",
                "Subrealm": "vira", # extra field
                "Kingdom": "virae",
                "Subkingdom": "virites", # extra field
                "Phylum": "viricota",
                "Subphylum": "viricotina", # extra field
                "Class": "viricetes",
                "Subclass": "viricetidae", # extra field
                "Order": "virales",
                "Suborder": "virineae", # extra field
                "Family": "viridae",
                "Subfamily": "virinae",
                "Genus": "virus",
                "Subgenus": "virus" # again! Extra field. Not used below beacause it's indistinguishable from "Genus".
                } # "Species" has no suffix (see above).
    

    # for every row in csv_edit (that means: for every row in the vConTACT2 csv output)
    for index, row in csv_edit.iterrows():
        
        # find the corresponding GenBank accession in the metas table:
        matches = metas[metas['Accession'] == row.Accession]
        if len(matches) == 0:
            # so this is a vOTU (that means: a viral scaffold in need to be classified)
            pass # go on with the next row!
        elif (len(matches) != 0 and len(matches) != 1):
            # there should be only 1 exact match.
            consoleout('error', "More than 1 match in filler() function!") # stop the program.
        else:
            # get first (and only) match:
            match = matches.iloc[0]

            # fill csv_edit with missing infos (Species, Genus, Family, Order, ...) from match (metas)
            csv_edit.at[index,'Host'] = match.Host

            
            csv_edit.at[index,'BaltimoreGroup'] = match.BaltimoreGroup
            csv_edit.at[index,'Realm'] = match.Realm
            # Subreal # -
            # Kingdom # -
            # Subkingdom # -
            # Phylum # -
            # Subphylum # -
            # Class # -
            # Subclass # -
            # Order # - 
            # Suborder # - 
            csv_edit.at[index,'Family'] = match.Family 
            csv_edit.at[index,'Subfamily'] = match.Subfamily 
            csv_edit.at[index,'Genus'] = match.Genus
            csv_edit.at[index,'Species'] = match.Description # match.Description contains only the species name.
                        
            
            # Below we try to catch the missing isolated fileds from metas (Kindom, Phylum, Class, ...)
            # from match.Classification, that is a big string like containing the full known linalogy, for example:
            # "Pseudomonas virus phiCTX Citexvirus Peduovirinae Myoviridae Caudovirales Caudoviricetes Uroviricota Heunggongvirae Duplodnaviria Viruses"
            
            # Anyway, match.Classification (the linalogy) is rarely complete. Some levels are missing here or there.
            # Fortunately, for Kindom, Phylum, Class etc there is a "suffix convention":
            # eg: _virae for Kingdom, _viricota for Phylum, _viricetes for Class. E
            # Please note that very suffix was before saved in lvlsuffix.
            
            # Remove the last level (species) from linealogy, then split by ' ':
            linealogy = match.Classification.replace(match.Description + " ","").split(" ")

            def scroll_linealogy(lin, suff): # function to find term with correspondent suffix in linealogy (list).
                cnt = 0
                got = "n.a." # Maybe it could be better to replace "n.a." with somthing more specific.
                for x in lin:
                    if x.endswith(suff):
                        cnt += 1
                        if cnt > 1: # There SHOULD be only one match for every suffix:
                            consoleout('warning', "More than 1 term in linealogy with suffix " +suff+ "! Taking the last.") 
                        got = x
                return got
            
            csv_edit.at[index,'Subrealm'] = scroll_linealogy(linealogy, lvlsuffix['Subrealm'])
            csv_edit.at[index,'Kingdom'] = scroll_linealogy(linealogy, lvlsuffix['Kingdom'])
            csv_edit.at[index,'Subkingdom'] = scroll_linealogy(linealogy, lvlsuffix['Subkingdom'])
            csv_edit.at[index,'Phylum'] = scroll_linealogy(linealogy, lvlsuffix['Phylum'])
            csv_edit.at[index,'Subphylum'] = scroll_linealogy(linealogy, lvlsuffix['Subphylum'])
            csv_edit.at[index,'Class'] = scroll_linealogy(linealogy, lvlsuffix['Class'])
            csv_edit.at[index,'Subclass'] = scroll_linealogy(linealogy, lvlsuffix['Subclass'])
            csv_edit.at[index,'Order'] = scroll_linealogy(linealogy, lvlsuffix['Order'])
            csv_edit.at[index,'Suborder'] = scroll_linealogy(linealogy, lvlsuffix['Suborder'])

            # Anyway store the full linealogy for following eventual uses.
            csv_edit.at[index,'FullClass'] = match.Classification

            # Note: after running the script, seems like metas 1Nov2021_data_excluding_refseq.tsv doesn not contain
            # Subreal, Subkingdom, Subphylum, Subclass and Suborder informations in metas.Classification.

    # Saving csv_edit as Excel:
    csv_edit.to_excel(output_path + 'csv_edit_' + string_suffix + '.xlsx' )

    # Saving csv_edit object as binary file:
    with open(output_path + 'csv_edit_' + string_suffix + '.bin', 'wb') as csv_edit_file: # wb = write binary.
        pickle.dump(csv_edit, csv_edit_file)

    consoleout('okay', 'Finished to flush the INPHARED taxonomy.' )
    return csv_edit



if __name__ == "__main__":

    # This is how to test the program:
    """
    python  graphanalyzer.py \
    --graph ./testinput/c1.ntw \
    --csv   ./testinput/genome_by_genome_overview.csv \
    --metas ./testinput/1Nov2021_data_excluding_refseq.tsv \
    --output      ./testoutput/ \
    --suffix      assemblerX \
    --preprocessed ./testoutput/csv_edit_assemblerX.bin 
    """
    # get the parameters from argparser:
    parameters = parser.parse_args()

    # check the presence of passed files:
    try:
        graph_table = open(parameters.graph, 'r') 
    except RuntimeError:
        consoleout('error', "Can't find the --graph passed as %s." % parameters.graph)
    try:
        csv_table = open(parameters.csv, 'r') 
    except RuntimeError:
        consoleout('error', "Can't find the --csv passed as %s." % parameters.csv)
    try:
        metas_table = open(parameters.metas, 'r') 
    except RuntimeError:
        consoleout('error', "Can't find the --metas passed as %s." % parameters.metas)
    consoleout('okay', 'Each required input file seems correctly loaded.')
    # Now we want to check if the goodness of the filepath provided. 
    if (os.path.isdir(parameters.output) == False): # if it's not a folder:
        if (os.path.isfile(parameters.output) == True): # if it's file:
            consoleout('error', "Seems that --output passed as %s is a file and not a folder." % parameters.output)
        else: # it's a folder, but not yet created. 
            consoleout('error', "Seems that --output passed as %s doesn't exist." % parameters.output) 

    # Now we try to load the graph_table into a NetoworkX object. From many types, we choose Graph:
    """
    from: https://networkx.github.io/documentation/stable/reference/classes/index.html#which-graph-class-should-i-use
    Networkx Class      Type            Self-loops allowed      Parallel edges allowed
    Graph               undirected      Yes                     No
    DiGraph             directed        Yes                     No
    MultiGraph          undirected      Yes                     Yes
    MultiDiGraph        directed        Yes                     Yes
    """

    # Load the graph_table into a networkx.Graph(). 
    graph = net.read_edgelist(graph_table, nodetype=str, data=(('weight',float),), create_using=net.Graph())
    graph_table.close()

    # Calculate the arrows' weight distribution:
    graph_table = open(parameters.graph, 'r')
    arrows = pnd.read_csv(graph_table, header = None, sep=' ')
    global_weights = list(arrows[2])
    """
    print("global_weights len: ", len(global_weights))
    print("global_weights MAX: ", max(global_weights))
    print("global_weights min: ", min(global_weights))
    print("global_weights mean: ", stats.mean(global_weights))
    print("global_weights mode: ", stats.mode(global_weights))
    print("global_weights median: ", stats.median(global_weights))
    """
    
    # 1st PART:
    # Here we want to fill the vConTACT2 csv output with the metas (taxonomy, etc) provided by INPHARED:
    # Developer could already provide the csv_edit table, just for the sake of speeding up the execution time:
    if parameters.preprocessed != None: 
        with open(parameters.preprocessed, 'rb') as csv_edit_file: # rb = read binary file.
            csv_edit = pickle.load(csv_edit_file)
            consoleout("okay", "Skipping the first part of the script since a preprocessed genome_by_genome_overview.csv was provided.")
    else: # continue with the normal flow of the script.

        # load csv and metas as pandas dataframe:
        csv = pnd.read_csv(csv_table, header = 0)
        metas = pnd.read_csv(metas_table, header = 0, sep='\t')

        # These should be all our vOTUs (just a list of their names):
        votus = csv[csv['Genome'].str.contains('vOTU')]['Genome'].tolist()
        # These should be all vOTUs contained in the graph:
        votus_ingraph = [votu for votu in votus if graph.has_node(votu) == True ]
        # Here we want ot check how many vOTU are there in 'csv'.
        # This should be equivalent to 'grep -Eo "vOTU_.*?," genome_by_genome_overview.csv | sort | uniq | wc -l'
        consoleout("okay", "There are " + str(len(votus)) + " vOTU in input.")
        # Here we want ot check how many vOTU are there in 'metas'.
        # This should be equivalent to 'grep -Eo "vOTU_.*? " c1.ntw | sort | uniq | wc -l'
        consoleout("okay", str(len(votus_ingraph)) + " vOTU are contained in the graph.")

        # fill the vConTACT2 csv output with the taxonomy provided by INPHARED:
        csv_edit = fillWithMetas(csv, metas, parameters.output, parameters.suffix)

    
    # 2nd PART:
    # Run the main algorithm. For every viral scaffold, get the most probable taxonomy: 
    image_table, df_results = clusterExtractor(graph, csv_edit, parameters.output, parameters.suffix)
    consoleout("okay", "Finished to assign taxonomy to vOTUs.")
    
    
    # 3rd PART:
    # generate a general plot with all reference genomes and viral scaffolds:
    image_graph = plotCreatorGraphvizHoloviews(graph, csv_edit, df_results, parameters.output, parameters.suffix, max(global_weights)) 


    # 4th PART:
    # MAKING OF THE PANEL APP (https://holoviz.org/tutorial/Building_Panels.html)
    # pane: view of an external object (text, image, plot, etc.) by wrapping it.
    # panel: lays out multiple components in a row, column, or grid.
    # widget: provides input controls to add interactive features to the panel.
    title_pane = pnl.panel('<h1 style="text-align: center;"> ' + parameters.suffix + ' graph panel </h1>', width=800)
    graph_pane = pnl.panel(image_graph) # returns a HoloViews(Overlay)
    table_pane = pnl.panel(image_table) # returns a HoloViews(Table)
    columnar_panel = pnl.Column(title_pane, graph_pane, table_pane)
    columnar_panel.save(parameters.output + 'panel_graph_' + parameters.suffix + '.html', resources=INLINE)
    consoleout("okay", "Finished to generate the panel app.")