# Coded by Gioele Lazzari (gioele.lazza@studenti.univr.it)
software = "graph_analyzer.py"
version = "0.8.3"

import sys, os, io
import argparse
import logging
from operator import itemgetter
from copy import deepcopy


#############
# ARGPARSER #
#############

parser = argparse.ArgumentParser(
    description = "This script is designed to extract the most probable taxonomic "
                  "assignation from the genome_by_genome_overview.csv file generated "
                  "by vConTACT2. In fact, on bitbucket.org/MAVERICLab/vcontact2/ "
                  "at date 21-09-2020, it's possible to read what follows: \n"
                  "One important note is that the taxonomic information is not included "
                  "for user sequences. This means that each user will need to find their "
                  "genome(s) of interest and check to see if reference genomes are located "
                  "in the same VC. If the user genome is within the same VC subcluster as a "
                  "reference genome, then there's a very high probability that the user "
                  "genome is part of the same genus. If the user genome is in the same VC "
                  "but not the same subcluster as a reference, then it's highly likely the two "
                  "genomes are related at roughly genus-subfamily level. If there are no "
                  "reference genomes in the same VC or VC subcluster, then it's likely that "
                  "they are not related at the genus level at all. That said, it is possible "
                  "they could be related at a higher taxonomic level (subfamily, family, order).",
    formatter_class = argparse.RawTextHelpFormatter)

options = parser.add_argument_group("Options")

options.add_argument('-v', '--version', action='version', version= software + " v" + version)

options.add_argument('-g', '--input-graph', dest='input_graph', metavar='FILENAME', default='c1.ntw')
options.add_argument('-c', '--input-csv', dest='input_csv', metavar='FILENAME', default='genome_by_genome_overview.csv')
options.add_argument('-o', '--output', dest='output_path', metavar='PATH', default='./')
options.add_argument('-p', '--suffix', dest='string_suffix', metavar='STRING', default='sample')


######################
# EXTERNAL LIBRARIES #
######################

def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)

try:
    import pandas as pnd 
except ImportError:
    error("The pandas library was not found.")

try:
    import networkx as net # net.spring_layout() could require also scipy
except ImportError:
    error("The networkx or scipy library was not found.")

try: 
    from networkx.drawing.nx_agraph import graphviz_layout
    # Mac: conda install graphviz==2.42.3=h055b950_2 -c conda-forge
    # Mac: conda install pygraphviz==1.6=py37hd4be98e_1 -c conda-forge
except ImportError:
    error("The pygraphviz library was not found.") 

try:
    import hvplot.networkx as hvnx
    import hvplot.pandas # hvPlot dynamically adds the Pandas .hvplot() method
    import hvplot # for the save() method
    from bokeh.resources import INLINE # for viewing html pages also without an internet connection
except ImportError:
    error("The hvplot library was not found.")

try:
    import panel as pnl
    pnl.extension() # before displaying anything with Panel it is always necessary to load the Panel extension
except ImportError:
    error("The panel library was not found.")


##################
# CORE FUNCTIONS #
##################

def clusterExtractor(g, df, output_path, string_suffix):

    def insertReference(df, results, scaffold, best_reference, weight_heavier_reference, level_reference, status):

        # extract taxonomy
        order = "-"
        family = "-"
        genus = "-"
        if level_reference != "F" and level_reference != "U":
            order = df.loc[df['Genome'] == best_reference, 'Order'].values[0]
            family = df.loc[df['Genome'] == best_reference, 'Family'].values[0]
            genus = df.loc[df['Genome'] == best_reference, 'Genus'].values[0]
            
        
        results = results.append({
            'Scaffold': scaffold,
            'status': status,
            'best_reference': best_reference.replace("~", " "),
            'level_reference': level_reference,
            'weight_reference': weight_heavier_reference,
            'Order': order,
            'Family': family,
            'Genus': genus,
            }, ignore_index = True) 

        return results

    ####################
    # PREPARATORY PART #
    ####################

    # create the .log file
    logging.basicConfig(filename = output_path + 'graph_analyzer_' + string_suffix + '.log', filemode='w', level = logging.INFO, # level sets the threshold
                        format = '%(asctime)s %(levelname)s: %(message)s',
                        datefmt = '%H:%M:%S') 
    logging.info('Processing of vConTACT2 output is started.\n')

    # less problematic when called as variable
    df = df.rename(columns={"VC Subcluster": "VCSubcluster"})
    df = df.rename(columns={"VC Status": "VCStatus"})

    # sobstitute NA with '' to avoid future runtime errors
    df = df.fillna('')

    # make a copy of the original df
    df_copy = df.copy(deep = True)
    # add a column that is a copy of "Genome": it will be updated at every iteration of the algorithm
    df_copy['Genome_editable'] = df['Genome']
    # make a copy of the original graph. This will cntain an editable flag for each scaffold node
    g_copy = deepcopy(g)
    for node in g_copy.nodes: # create an "assignment" attribute for all nodes
        net.set_node_attributes(g_copy, {node: node}, "assignment") # store a copy of node's name

    # prepare results dataframe
    results = pnd.DataFrame(data = {
        'Scaffold': [],
        'status': [],
        'best_reference': [],
        'level_reference': [],
        'weight_reference': [],
        'Order': [],
        'Family': [],
        'Genus': []})

    ###################
    # PROCESSING PART #
    ###################

    # extract scaffolds from reference genomes
    scaffolds_total = df_copy[df_copy['Genome_editable'].str.contains('NODE')]
    scaffolds_ingraph = scaffolds_total.copy(deep = True)
    for row in scaffolds_ingraph.itertuples():
        if g_copy.has_node(row.Genome) == False : # Note: not all scaffolds are present in the graph!
            # remove row if scaffold is not in the graph
            scaffolds_ingraph = scaffolds_ingraph.drop(scaffolds_ingraph[scaffolds_ingraph.Genome == row.Genome].index)
            # append in results as level_reference='U'(UNCLSUTERED)
            results = insertReference(df=df, results=results, scaffold=row.Genome, best_reference="-", weight_heavier_reference=0.0, level_reference='U', status="-")
                    
    logging.info("Viral scaffolds in total: %i" % len(scaffolds_total))
    logging.info("Viral scaffolds in the graph: %i\n" % len(scaffolds_ingraph))



    # to count how many iterations the algorithm does
    counter_iterations = 0
    # to count how many new assignments at every iteration
    counter_new = -1 # -1 is just to start the algorithm
    while(counter_new != 0):
        counter_new = 0 # reset the counter at every iteration
        counter_iterations += 1
        logging.info("##############################################")
        logging.info("################ Iteration %i #################" % counter_iterations)
        logging.info("##############################################\n")

        
        # for each viral scaffold:
        for row in scaffolds_ingraph.itertuples():
            
            # skip already assigned scaffolds
            if not "NODE" in row.Genome_editable:
                continue
            
            # extract scaffold name and VC
            scaffold = row.Genome
            status = row.VCStatus
            logging.info("Scaffold: %s" % (scaffold))
            logging.info("VCStatus: %s" % (status))
              
            # COMPUTE NEIGHBORS
            neighbors_list = list(g_copy.neighbors(scaffold))
            logging.info("Neighbors: %i" % len(neighbors_list))
            

            if (status == "Clustered" or status == "Outlier" or status == "Clustered/Singleton" or "Overlap" in status):
              
                vc = "" # scaffold's viral cluster or subcluster
                sameclustered_list = [] # other genomes/scaffolds in the same cluster
                
                if status == "Clustered":
                    # EXTRACT SAMECLUSTERED
                    vc = row.VCSubcluster
                    # extract all scaffolds clustered in the same subcluster (hereafter: the sameclustered)
                    sameclustered = df_copy[df_copy['VCSubcluster'].str.contains(vc)]
                    # remove current scaffold's row
                    sameclustered = sameclustered.drop(sameclustered[sameclustered.Genome == scaffold].index)
                    sameclustered_list = sameclustered["Genome"].tolist()
                
                elif status == "Outlier": # not clustered but connected
                    vc = ""
                    sameclustered_list = []

                elif status == "Clustered/Singleton": # they have a subcluster of their own
                    # extract all the possible subclusteres within the cluster
                    vc = "VC_" + row.VCSubcluster.split("_")[1]
                    sameclustered = df_copy[df_copy['VCSubcluster'].str.contains(vc)]
                    # remove current scaffold's row
                    sameclustered = sameclustered.drop(sameclustered[sameclustered.Genome == scaffold].index)
                    sameclustered_list = sameclustered["Genome"].tolist()
                    

                elif "Overlap" in status: # "Overlap" belong to n cluster
                    # extract the n cluster and consider every of thier subcluster
                    vc_list = status.replace("Overlap (", "").replace(")", "").split("/")
                    sameclustered = pnd.DataFrame() # empty df
                    for vc in vc_list: # extract the rows and glue them the the previous
                        sameclustered = sameclustered.append(df_copy[df_copy['VCSubcluster'].str.contains(vc)])
                    # remove current scaffold's row
                    sameclustered = sameclustered.drop(sameclustered[sameclustered.Genome == scaffold].index)
                    sameclustered_list = sameclustered["Genome"].tolist()


                logging.info("Sameclustered: %i" % len(sameclustered_list))
                # check if all sameclustered are contained in the neighbours
                logging.info("All connected?: %s" % set(sameclustered_list).issubset(neighbors_list))

                
                # GET SAMECLUSTERED-CONNECTED LIST ORDERED BY WEIGHT
                # associate every CONNECTED sameclustered with its weight
                connected_list = [] # this will be a list of dict
                for node in sameclustered_list:
                    try: 
                        w = g_copy.get_edge_data(scaffold, node)["weight"]
                        assignment = g_copy.nodes[node]['assignment']
                        connected_list.append({"scaffold": node, "weight": w, "assignment": assignment})
                    except: # not every pair of nodes in the same cluster are directly connected in the graph
                        continue
                # ordinate the list by weight (itemgetter is from the operator library)
                connected_list = sorted(connected_list, key=itemgetter('weight'), reverse=True) # reverse for descending order


                # EXTRACT BEST_BET (= the haviest connected-sameclustered)
                # note that best_bet could be another scaffold!
                try: 
                    best_bet = connected_list[0]["scaffold"] # to store the node connected with the heavier edge
                    weight_heavier_edge = connected_list[0]["weight"] # to store the heavier edge
                except: # not always there is a best_bet: for example "Outlier" are not clustered
                    best_bet = "-"
                    weight_heavier_edge = 0.0 
                logging.info("best_bet: %s" % best_bet)
                logging.info("weight_heavier_edge: %f" % weight_heavier_edge)

                
                # BEST_REFERENCE SEARCH
                best_reference = "-" # store the reference genome connected with the heavier edge
                weight_heavier_reference = 0.0 # store the heavier egde that connect to a reference genome
                level_reference = "-" # keep track of the cardinality (order of the list)
                # iterate the list until the first reference genome is reached
                for i in range(len(connected_list)):
                    # if it's connected a true reference genome or a scaffold assigned in the last iteration:
                    if  connected_list[i]["scaffold"].find("NODE") == -1 or connected_list[i]["assignment"].find("NODE") == -1:
                        
                        if connected_list[i]["scaffold"].find("NODE") == -1:
                            best_reference = connected_list[i]["scaffold"]
                            weight_heavier_reference = connected_list[i]["weight"]
                            level_reference = str(counter_iterations) + "C" + str(i + 1) # first level is 1, not 0

                        elif connected_list[i]["assignment"].find("NODE") == -1 and counter_iterations == 1:
                            continue # scaffold assigned are not considerated in the first iteration. This will skip the "break"

                        elif connected_list[i]["assignment"].find("NODE") == -1:
                            best_reference = connected_list[i]["assignment"]
                            weight_heavier_reference = connected_list[i]["weight"]
                            level_reference = str(counter_iterations) + "C" + str(i + 1) # first level is 1, not 0
                        
                        # break the for-loop, because the first reference genome was reached
                        break
                        

                # if a reference genome was found within the connected-sameclustered:            
                if best_reference != "-": 
                    counter_new += 1

                    logging.info("VCSubcluster: %s" % (vc))
                    logging.info("best_reference: %s" % best_reference)
                    logging.info("weight_heavier_reference: %f" % weight_heavier_reference) 
                    logging.info("level_reference: %s" % level_reference)          

                    # update results
                    results = insertReference(df, results, scaffold, best_reference, weight_heavier_reference, level_reference, status)
                    # upadate starting table and graph
                    scaffolds_ingraph.loc[scaffolds_ingraph["Genome"] == scaffold, ["Genome_editable"]] = best_reference
                    g_copy.nodes[scaffold]['assignment'] = best_reference
                
                
                else: # RETRY ITERATING THE WHOLE NEGHBORS 
                    
                    # associate every neighbors with its weight
                    connected_list = [] # this will be a list of dict
                    for node in neighbors_list:
                        w = g_copy.get_edge_data(scaffold, node)["weight"]
                        assignment = g_copy.nodes[node]['assignment']
                        connected_list.append({"scaffold": node, "weight": w, "assignment": assignment})
                    # ordinate the list by weight (itemgetter is from the operator library)
                    connected_list = sorted(connected_list, key=itemgetter('weight'), reverse=True) # reverse for descending order
                    
                    
                    # BEST_REFERENCE SEARCH
                    # iterate the list until the first reference genome is reached
                    for i in range(len(connected_list)):
                        if  connected_list[i]["scaffold"].find("NODE") == -1 or connected_list[i]["assignment"].find("NODE") == -1:

                            if connected_list[i]["scaffold"].find("NODE") == -1:
                                best_reference = connected_list[i]["scaffold"]
                                weight_heavier_reference = connected_list[i]["weight"]
                                level_reference = str(counter_iterations) + "N" + str(i + 1) # first level is 1, not 0

                            elif connected_list[i]["assignment"].find("NODE") == -1 and counter_iterations == 1:
                                continue # scaffold assigned are not considerated in the first iteration. This will skip the "break"
                            
                            elif connected_list[i]["assignment"].find("NODE") == -1:
                                best_reference = connected_list[i]["assignment"]
                                weight_heavier_reference = connected_list[i]["weight"]
                                level_reference = str(counter_iterations) + "N" + str(i + 1) # first level is 1, not 0
                            
                            # break the for-loop, because the first reference genome was reached
                            break


                    # if a reference genome was found within the whole neighbors:                                 
                    if best_reference != "-": 
                        counter_new += 1
                        
                        logging.info("VCSubcluster: %s" % (vc))
                        logging.info("best_reference: %s" % best_reference)
                        logging.info("weight_heavier_reference: %f" % weight_heavier_reference) 
                        logging.info("level_reference: %s" % level_reference)          

                        # update results
                        results = insertReference(df, results, scaffold, best_reference, weight_heavier_reference, level_reference, status)
                        # upadate starting table and graph
                        scaffolds_ingraph.loc[scaffolds_ingraph["Genome"] == scaffold, ["Genome_editable"]] = best_reference
                        g_copy.nodes[scaffold]['assignment'] = best_reference


                    else: # this means: NO reference in sameclustered AND NO reference in neighbors
                        logging.warning("Iteration %i: NO reference in connected sameclustered AND NO reference in neighbors." % counter_iterations)
               
            else: # "Singleton" never present in the graph
                logging.error("Something else happened!")   

            logging.info("##############################################\n")

    # add the remaining scaffolds (level_reference='F': present in the graph but not assigned)
    for row in scaffolds_ingraph.itertuples():
        if "NODE" in row.Genome_editable: 
            results = insertReference(df=df, results=results, scaffold=row.Genome_editable, best_reference="-", weight_heavier_reference=0.0, level_reference='F', status="-")

    logging.info('Processing of vConTACT2 output is ended.') 
    
    # order results by scaffold name
    results = results.sort_values(by=['Scaffold'])
    # convert the dataframe to a html table
    image_table = results.hvplot.table(width = 800, height = 200)
    # save results as a csv table
    results.to_csv(output_path + "results_vcontact2_" + string_suffix + ".csv")
    # save results as a MultiQC custom table
    content = """# plot_type: 'table'
# section_name: 'vConTACT2 taxonomy table'
# description: 'Taxonomy table: automatic processing of vConTACT2 outputs `c1.ntw` and `genome_by_genome_overview.csv`.'
<-- REPLACE -->
"""
    string_eater = io.StringIO()
    datas = results.to_csv( path_or_buf = string_eater, sep='\t', index=False)
    content = content.replace("<-- REPLACE -->", string_eater.getvalue())
    file_report = open("custom_taxonomy_table_mqc.txt", "w")
    file_report.write(content)
    file_report.close()

    return image_table, results



def plotCreatorGraphvizHoloviews(g, df, output_path, string_suffix, results):

    ####################
    # PREPARATORY PART #
    ####################

    # drop from results every scaffold not assigned (FREE or UNCLUSTERED)
    results_copy = deepcopy(results)
    for row in results.itertuples():
        if row.level_reference == "F" or row.level_reference == "U": 
            results_copy = results_copy.drop(results_copy[results_copy.Scaffold == row.Scaffold].index)
            

    ###################
    # PROCESSING PART #
    ###################

    # extract viral scaffolds from reference genomes (only the scaffolds present in the graph)
    scaffolds, present_nodes = df[df['Genome'].str.contains('NODE')] , []
    for row in scaffolds.itertuples():
        if g.has_node(row.Genome): # Note: not all scaffolds are clustered!
            present_nodes.append(row.Genome)
    
    # extract scaffold that have received a taxanomy
    assigned_nodes = results_copy['Scaffold'].tolist()
    # extract all nodes connected with "presents" nodes
    connected_nodes = []
    for element in present_nodes:
        connected_nodes = connected_nodes  + list(net.node_connected_component(g, element))
    connected_nodes = list(set(connected_nodes)) # remove duplicates if present
    connected_g = g.subgraph(connected_nodes) # create a sugraph for speeding up subsequent calcs
    # exclude assigned_nodes from the present_nodes
    present_nodes = list(set(present_nodes)-set(assigned_nodes))
    # exclude scaffold from the reference genomes
    remaining_nodes = list(set(connected_nodes)-set(present_nodes)-set(assigned_nodes))
    

    # LABEL CREATION
    labels = results_copy['best_reference'].tolist()
    # create dictionary for labels
    labdict = {}
    for i in range(len(assigned_nodes)):
        labdict[assigned_nodes[i]] = labels[i]


    # calculate position of  nodes
    # K: roughly corresponds to an ideal edge length (in inches). len can be used to override this value for adjacent nodes.
    # repulsiveforce: values larger than 1 tend to reduce the warping effect at the expense of less clustering.
    # overlap: determines if and how node overlaps should be removed. If "true" , overlaps are retained. 
    pos = graphviz_layout(connected_g, prog="sfdp",  args='-Goverlap=true')
    # note: sfdp requires graphviz built with gts, but this type of build is rarely available in conda
    # so sfdp will print "Error: remove_overlap: Graphviz not built with triangulation library"
    # right builds are graphviz==2.42.3=h055b950_2 and pygraphviz==1.6=py37hd4be98e_1 from conda-forge


    # PLOTTING as described in https://hvplot.holoviz.org/user_guide/NetworkX.html 
    # colors available at https://docs.bokeh.org/en/latest/docs/reference/colors.html 
    # in holoviews and in hvplot, the * operator superimposes the different layers of the plot
    # first of all, draw all the edges with the same style
    image = hvnx.draw_networkx_edges(connected_g, pos=pos, edge_width = 0.1, alpha = 0.3, width = 800, height = 500)
    # then draw the reference genomes 
    image = image * hvnx.draw_networkx_nodes(g.subgraph(remaining_nodes), pos=pos, node_color="aquamarine", node_size=200, alpha = 0.3, linewidths =0.0) 
    # draw scaffold not assigned
    image = image * hvnx.draw_networkx_nodes(g.subgraph(present_nodes), pos=pos, node_color="blue", node_size=200, alpha = 0.5, linewidths =0.0) 
    # draw assigned scaffolds with labels
    image = image * hvnx.draw_networkx_nodes(g.subgraph(assigned_nodes), pos=pos, node_color="red", node_size=200, alpha = 0.5,  linewidths =0.0, labels= labdict, font_size = "8pt")
    
    return image



if __name__ == "__main__":

    terminal_mode = True
    if not terminal_mode:
        sys.argv = ['graph_analyzer.py', 
                '--input-graph', './inputs/c1.ntw',
                '--input-csv', './inputs/genome_by_genome_overview.csv', 
                '--output', './outputs/',
                '--suffix', 'crosta']
    parameters = parser.parse_args()


    with open(parameters.input_graph, 'r') as graph_table ,  open(parameters.input_csv, 'r') as csv_table:
        
        """
        from: https://networkx.github.io/documentation/stable/reference/classes/index.html#which-graph-class-should-i-use
        Networkx Class      Type            Self-loops allowed      Parallel edges allowed
        Graph               undirected      Yes                     No
        DiGraph             directed        Yes                     No
        MultiGraph          undirected      Yes                     Yes
        MultiDiGraph        directed        Yes                     Yes
        """

        # this is the total graph, usually opened with Cytoscape
        g = net.read_edgelist(graph_table, nodetype=str, data=(('weight',float),), create_using=net.Graph())

        # this is the table that specifies the clusters division
        df = pnd.read_csv(csv_table, header = 0)

        # for every viral scaffold, get the most probable taxonomy 
        image_table, results = clusterExtractor(g, df, parameters.output_path, parameters.string_suffix)
        
        # generate a general plot with all reference genomes and viral scaffolds
        image_graph = plotCreatorGraphvizHoloviews(g, df, parameters.output_path, parameters.string_suffix, results) 

        # MAKING OF THE PANEL APP (https://holoviz.org/tutorial/Building_Panels.html)
        # pane: view of an external object (text, image, plot, etc.) by wrapping it
        # panel: lays out multiple components in a row, column, or grid
        # widget: provides input controls to add interactive features to the panel
        title_pane = pnl.panel('<h1 style="text-align: center;"> ' + parameters.string_suffix + ' graph panel </h1>', width=800)
        graph_pane = pnl.panel(image_graph) # returns a HoloViews(Overlay)
        table_pane = pnl.panel(image_table) # returns a HoloViews(Table)
        columnar_panel = pnl.Column(title_pane, graph_pane, table_pane)
        columnar_panel.save(parameters.output_path + 'panel_graph_' + parameters.string_suffix + '.html', resources=INLINE)
        
    