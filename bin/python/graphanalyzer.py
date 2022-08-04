#!/usr/bin/env python3
# Coded by Gioele Lazzari (gioele.lazzari@univr.it)
software = "graphanalyzer.py"
version = "1.4" 


# import system libraries:
import sys, os, io, argparse, logging, textwrap, time, datetime
from operator import itemgetter
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor, as_completed
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
        sys.stdout.write(msg)
        sys.stdout.flush()
    elif lvl == 'okay': 
        msg = colors['OKGREEN'] + 'OKAY: ' + colors['ENDC'] + msg + '\n'
        sys.stdout.write(msg)
        sys.stdout.flush()
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
parser.add_argument('-p', '--prefix', 
                    metavar='STRING', 
                    default='vOTU',  
                    type=str,
                    help='prefix of the header of each conting representing a vOTU')
parser.add_argument('-s', '--suffix', 
                    metavar='STRING', 
                    default='assemblerX',  
                    type=str,
                    help='suffix to append to every file produced in the output directory')
parser.add_argument('-t', '--threads',
                    metavar='INT', 
                    default=4,  
                    type=int,
                    help='how many threads to use for the generation of the interactive subgraphs')


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






 
def fillWithMetas(csv , metas):

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
    # "Genus", "Family", and "Order" are commented as they were already created in the original "csv".
    # "Genus" and "Family" will anyway be filled with info from "metas", matching the accession number. 

    # Taxonomic levels signed with *, are not present in "metas" tables. We will try to fill
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
    # Order     # *
    csv_edit["Suborder"] = None # *
    # Family      # to be filled with info from metas, matching the accession number.
    csv_edit["Subfamily"] = None  # to be filled with info from metas, matching the accession number.
    # Genus       # to be filled with info from metas, matching the accession number.
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
    

    # for every row in csv_edit (that means: for every row in the vConTACT2 "csv" output)
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
                            consoleout('warning', f'More than 1 term in linealogy with suffix "{suff}"! Taking the last.') 
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

    return csv_edit



def clusterExtractor(graph, csv_edit, output_path, string_suffix, prefix):

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
        # "F": present in the graph but not assigned (completely isolated).
        # "A": present in the graph but not assigned (indireclty connected to references).
        if level != "F" and level != "G" and level != "A":
            
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


    ###################
    # PROCESSING PART #
    ###################

    # extract scaffolds from reference genomes:
    scaffolds_total = csv_edit[csv_edit['Genome'].str.startswith(prefix)]

    # extract scaffolds contained in the graph (as pnd dataframe):
    scaffolds_ingraph = scaffolds_total.copy(deep = True)
    for index, row in scaffolds_total.iterrows():
        if graph.has_node(row.Genome) == False : 
            
            # remove row from scaffolds_ingraph if scaffold is not in the graph:
            scaffolds_ingraph = scaffolds_ingraph.drop(index)
            
            # In this case "status" should be always "Singleton" (that is: never present in the graph).
            status = row.VCStatus # SINGLETON

            # append dropped row in results as level="G" (meaning: NOT IN GRAPH):
            results = insertReference(
                csv_edit=csv_edit, 
                results=results, 
                scaffold=row.Genome, 
                closer="n.a.", 
                weight="n.a.", 
                level="G", 
                status=status, 
                vc="n.a."
            )

    # log first stats:                
    logging.info("Viral scaffolds in total: %i" % len(scaffolds_total))
    logging.info("Viral scaffolds in the graph: %i\n" % len(scaffolds_ingraph))

   
    # for each viral scaffold in the graph:
    for index, row in scaffolds_ingraph.iterrows():

        
        # extract scaffold name and VCStatus:
        scaffold = row.Genome
        status = row.VCStatus
        logging.info("Scaffold: %s" % (scaffold))
        logging.info("VCStatus: %s" % (status))
            
        # COMPUTE NEIGHBORS (that are: the directly connected nodes):
        neighbors_list = list(graph.neighbors(scaffold))
        logging.info("# Neighbors: %i" % len(neighbors_list))
        

        if (status == "Clustered" or status == "Outlier" or status == "Clustered/Singleton" or "Overlap" in status):
            
            vc = "" # scaffold's viral cluster or subcluster
            sameclustered_list = [] # other genomes/scaffolds in the same cluster/subcluster
            
            # EXTRACT SAMECLUSTERED LIST:

            if status == "Outlier": # not clustered but connected
                vc = "n.a."
                sameclustered_list = []

            elif status == "Clustered":
                vc = row.VCSubcluster
                # extract all scaffolds clustered in the same subcluster (hereafter: the "sameclustered")
                sameclustered = csv_edit[csv_edit['VCSubcluster'] == vc]
                sameclustered_list = sameclustered["Genome"].tolist()
                sameclustered_list.remove(scaffold) # remove current scaffold

            elif status == "Clustered/Singleton": # they have a subcluster of their own
                # extract all the possible subclustered within the cluster
                vc = "VC_" + row.VCSubcluster.split("_")[1] + "_"
                sameclustered = csv_edit[csv_edit['VCSubcluster'].str.startswith(vc)]
                sameclustered_list = sameclustered["Genome"].tolist()
                sameclustered_list.remove(scaffold) # remove current scaffold

            elif "Overlap" in status: # "Overlap" belong to n cluster
                # extract the n cluster and consider every of thier subcluster
                vc_list = status.replace("Overlap (", "").replace(")", "").split("/")
                sameclustered = pnd.DataFrame() # empty df
                for vc in vc_list: # extract the rows and glue them together with the others
                    sameclustered = sameclustered.append(csv_edit[csv_edit['VCSubcluster'].str.startswith(vc + "_")])
                    # handling of others Overlaps:
                    sameclustered = sameclustered.append(csv_edit[
                        (csv_edit['VCStatus'].str.contains('\(' + vc + '/')) | \
                        (csv_edit['VCStatus'].str.contains('/' + vc + '/')) | \
                        (csv_edit['VCStatus'].str.contains('/' + vc + '\)')) ])
                sameclustered_list = sameclustered["Genome"].tolist() 
                sameclustered_list.remove(scaffold) 

            else:
                consoleout("error", "Strange level found during the extraction of the cluster.")

            # log some stats
            logging.info("Status: %s" % status)
            logging.info("Number of sameclustered (sc): %i" % len(sameclustered_list))
            # check if all sameclustered are contained in the neighbours
            logging.info("Are all sc connected to scaffold?: %s" % set(sameclustered_list).issubset(neighbors_list))

            
            # GET SAMECLUSTERED-CONNECTED LIST ORDERED BY WEIGHT
            # associate every CONNECTED sameclustered with its weight
            connected_list = [] # this will be a list of DICT:
            # {"scaffold": ___, "weight": ___}
            for node in sameclustered_list:
                try: 
                    w = graph.get_edge_data(scaffold, node)["weight"]
                    connected_list.append({"scaffold": node, "weight": w})
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
                if  connected_list[i]["scaffold"].startswith(prefix) == False:
                    closer_ref = connected_list[i]["scaffold"]
                    closer_ref_w = connected_list[i]["weight"]
                    level = "C" + str(i + 1) # first level is 1, not 0
                    # break the for-loop, because the first reference genome was reached
                    break
                    
            # if a reference genome was found within the connected-sameclustered:            
            if closer_ref != "n.a.": 
                logging.info("closer_ref: %s" % closer_ref)
                logging.info("closer_ref_w: %f" % closer_ref_w) 
                logging.info("level: %s" % level)          

                # update results:
                if status == "Clustered/Singleton":
                    # 'row.VCSubcluster' instead of 'vc', beacuse 'vc' changed during the algorithm
                    results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, row.VCSubcluster)
                elif "Overlap" in status:
                    # 'n.a.' instead of 'vc', beacuse 'vc' changed during the algorithm
                    results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc="n.a.")
                else:
                    results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc)
                
            
            
            else: # RETRY ITERATING THE WHOLE NEGHBORS 
                
                # associate every neighbors with its weight
                connected_list = [] # this will be a list of DICT:
                # {"scaffold": ___, "weight": ___}
                for node in neighbors_list:
                    w = graph.get_edge_data(scaffold, node)["weight"]
                    connected_list.append({"scaffold": node, "weight": w})
                # ordinate the list by weight (itemgetter is from the operator library)
                connected_list = sorted(connected_list, key=itemgetter('weight'), reverse=True) # reverse for descending order
                
                # 'CLOSER_REF' SEARCH
                # iterate the list until the first reference genome is reached
                for i in range(len(connected_list)):
                    if  connected_list[i]["scaffold"].startswith(prefix) == False:
                        closer_ref = connected_list[i]["scaffold"]
                        closer_ref_w = connected_list[i]["weight"]
                        level = "N" + str(i + 1) # first level is 1, not 0
                        # break the for-loop, because the first reference genome was reached
                        break

                # if a reference genome was found within the whole neighbors:                                 
                if closer_ref != "n.a.":
                    logging.info("closer_ref: %s" % closer_ref)
                    logging.info("closer_ref_w: %f" % closer_ref_w) 
                    logging.info("level: %s" % level)          

                    # update results
                    if status == "Clustered/Singleton":
                        # 'row.VCSubcluster' instead of 'vc', beacuse 'vc' changed during the algorithm
                        results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, row.VCSubcluster)
                    elif "Overlap" in status:
                        # 'n.a.' instead of 'vc', beacuse 'vc' changed during the algorithm
                        results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc="n.a.")
                    else:
                        results = insertReference(csv_edit, results, scaffold, closer_ref, closer_ref_w, level, status, vc)


                else: # this means: NO reference in sameclustered AND NO reference in neighbors
                    logging.warning("NO reference in connected sameclustered AND NO reference in neighbors.")

                    currlev = "F" # F: completely isolated vOTUs; A: vOTUs distantly indirectly linked to references.
                    genomes_connected = list(net.node_connected_component(graph, row.Genome))
                    for genome in genomes_connected: # here "genome" is a string !
                        if genome.startswith(prefix)==False:
                            currlev = "A"
                            break
                        
                    # From what we've seen, status for "F" can be:
                    # Clustered, Outlier, Overlap. But NOT "Singleton" (that is: never present in the graph).
                    # Please note: also "Clustered/Singleton" is valid VCStatus, though never found in 'F' here till now.
                    status = row.VCStatus # should NOT be a "Singleton". 
                    vc = row.VCSubcluster
                    if status == "Outlier" or "Overlap" in status: 
                        vc = "n.a."
                    results = insertReference(csv_edit, results, scaffold=row.Genome, closer="n.a.", weight="n.a.", level=currlev, status=status, vc=vc)

            
        else: # "Singleton" never present in the graph
            logging.consoleout('error', "Found a strange VCStatus! Maybe it's a Singleton? Singletons should not be present inside the graph.")   

        logging.info("##############################################\n")

    logging.info('Processing of vConTACT2 output is ended.') 
    

    # order results by scaffold name
    results = results.sort_values(by=['Scaffold'])

    # Cut taxonomy at a level dependent of the "confidence" of the assignment:
    for index, row in results.iterrows():
        # We decided to stop at:
        # Genus for Subclusters;
        # Subfamily (if exists, otherwise Family) for Clusters;
        # Order for 'N' (Outliers).
        if row.Level.startswith('C') and row.Status == "Clustered":
            continue # the best case: keep everything
        elif row.Level.startswith('C') and (row.Status == "Clustered/Singleton" or "Overlap" in row.Status):
            results.at[index,'Genus'] = "O" # Please note: "O" stands for "omitted".
        elif row.Level.startswith('N') :
            results.at[index,'Genus'] = "O" 
            results.at[index,'Subfamily'] = "O"
            # res.at[index,'Family'] = "O"


    # Now we want to save results dataframe into 4 different formats:
    # Formats are: 1) csv file; 2) Excel file; 3) MultiQC custom table.
    # 1) save results as a csv file: 
    results.to_csv(output_path + "results_vcontact2_" + string_suffix + ".csv")
    # 2) save results as a csv file: 
    results.to_excel(output_path + "results_vcontact2_" + string_suffix + ".xlsx")
    # 3) Save results as a MultiQC custom table (return the object and save in 4th step):
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

    return results



def subgraph_generation(scaffold, desired_path):

    neigh = list(graph.neighbors(scaffold))
    # recreate a subgraph with: nieghbors + scaffold
    # we assume that neighbors() doesn't already include 'scaffold'.
    neigh.append(scaffold)
    sview_graph = graph.subgraph(neigh)
      

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
    sview_image = hvnx.draw_networkx(
        sview_graph, pos=sview_pos, 
        edgelist=[], nodelist=[],  
        width=1000, height=650
    )


    # update the attribute 'A6_Weight' for every node:
    for node in sview_graph: # here node is a string !
        attrs = {} # a dict.
        if node == scaffold: attrs = {node: {"A6_Weight": "origin"}}
        else: attrs = {node: {"A6_Weight": str(round(sview_graph[scaffold][node]["weight"],1))}}
        net.set_node_attributes(sview_graph, attrs)


    # Now we want to draw all the edges:
    # Drawing edges one-by-one is too slow and generates too heavy .html files.
    # Using 'cmap' and 'dim' is faster and lighter. For example 90 KB vs 1.3 MB for vOTU_1.
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
    <body><p style="text-align:center">Interactive plot generated with <strong>graphanalyzer.py</strong>. Please wait the loading.</p>
    <p style="text-align:center">User guide available at <a href="https://www.github.com/lazzarigioele/graphanalyzer/">github.com/lazzarigioele/graphanalyzer</a>.</p>
    <p style="text-align:center">Bugs can be reported to <a href= "mailto:gioele.lazzari@univr.it">gioele.lazzari@univr.it</a>.</p>
    """)
    # specify the version of GA:
    tags = tags.replace("graphanalyzer.py", f"graphanalyzer.py v{version}")
    file = open(desired_path + scaffold + '.html', "w")
    file.write(wholetext.replace("<body>", tags))
    file.close() # always close file streams!

    return scaffold



def subgraphCreator(graph, csv_edit, results, output_path, string_suffix, max_weight, prefix, nthreads):

    # less problematic when called as variable
    csv_edit = csv_edit.rename(columns={"VC Subcluster": "VCSubcluster"})
    csv_edit = csv_edit.rename(columns={"VC Status": "VCStatus"})

    # sobstitute NA with '' to avoid future runtime errors
    csv_edit = csv_edit.fillna('')


    # PART 1.
    # Below we add attributes to each node in the graph. Thanks to this, the interactive
    # visualization will be more informative (each attribute is shown when a node is hovered
    # with the mouse).
    attribs = {} # create empty dict:

    # for each row in csv_edit, preparate a col and a dict to rename nodes
    for index, row in csv_edit.iterrows():

        # not all references in genome_by_genome_overview.csv are contained in the graph!
        # for genomes/scaffolds not included in graph: ignore.
        if graph.has_node(row.Genome) == False:
            continue

        # don't want to rename scaffold at this point
        if row.Genome.startswith(prefix): # Discriminate accesions from vOTUs.

            # find the corresponding vOTU in the results table:
            matches = results[results['Scaffold'] == row.Genome]
            if len(matches) != 1:
                # there should be only 1 exact match
                consoleout("error", "We have len(matches) != 1 in subgraphCreator().")
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
    # Here we want to draw a neighbours-based plot for each vOTUs in graph.
    results_ingraph = deepcopy(results)    # not G
    for index, row in results.iterrows():
        if row.Level == "G":
            results_ingraph = results_ingraph.drop(index)
    # extract scaffold that are in-graph / that have received a taxanomy / that are confident:
    scaffolds_ingraph = results_ingraph['Scaffold'].tolist()

    # First of all create the subfolder for storing all the subplots.
    desired_path = output_path + 'single-views_' + string_suffix + '/'
    if (os.path.isdir(desired_path) == False): # if it's not been created yet:
        try:
            os.mkdir(desired_path) # make it!
        except:
            consoleout("error", "Can't create the output sub-folder '%s'. " % desired_path)
    
    # Suggested reading https://superfastpython.com/threadpoolexecutor-vs-processpoolexecutor/ 
    # on the difference between ThreadPoolExecutor and ProcessPoolExecutor. 
    # Remember: with precesses, functions are only picklable if they are defined at the top-level of a module.
    executor = ProcessPoolExecutor(nthreads) 
    futures = [] # list of tasks
    # For each 'scaffolds_ingraph', compute the nieghtbors:
    for scaffold in scaffolds_ingraph:
        futures.append(executor.submit(subgraph_generation, scaffold, desired_path))

    counter = 0
    # iterate over all submitted tasks and get results as they are available.
    # as_completed() is useful to manipulate results as they become available.
    # as_completed() return results in whatever order (as soon as they become available).
    for future in as_completed(futures):
        # get the result for the next completed task
        scaffold = future.result() # blocking call: wait for this task to completed
        counter += 1
        print(f"Completed subgraph {counter}/{len(scaffolds_ingraph)} ({scaffold})        ", end='\r')
    executor.shutdown() # blocking call: wait for all tasks to complete


    return None 



if __name__ == "__main__":

    # This is how to test the program:
    """
    python  graphanalyzer.py \
    --graph ./testinput/c1.ntw \
    --csv   ./testinput/genome_by_genome_overview.csv \
    --metas ./testinput/1Nov2021_data_excluding_refseq.tsv \
    --output      ./testoutput/ \
    --prefix      vOTU \
    --suffix      assemblerX \
    -t 4
    """

    print(f"Starting graphanalyzer.py v{version} on {datetime.datetime.now()}")
    
    # get the parameters from argparser:
    parameters = parser.parse_args()

    # check the presence of passed files:
    try: graph_table = open(parameters.graph, 'r') 
    except RuntimeError: consoleout('error', "Can't find the --graph passed as %s." % parameters.graph)
    try: csv_table = open(parameters.csv, 'r') 
    except RuntimeError: consoleout('error', "Can't find the --csv passed as %s." % parameters.csv)
    try: metas_table = open(parameters.metas, 'r') 
    except RuntimeError: consoleout('error', "Can't find the --metas passed as %s." % parameters.metas)
    consoleout('okay', 'Each required input file seems correctly loaded.')

    # Now we want to check the goodness of the filepath provided. 
    if (os.path.isdir(parameters.output) == False): # if it's not a folder:
        if (os.path.isfile(parameters.output) == True): # it's file:
            consoleout('error', "Seems that --output passed as %s is a file and not a folder." % parameters.output)
        else: # it's a folder, but not yet created. 
            consoleout('error', "Seems that --output passed as %s doesn't exist." % parameters.output) 
    consoleout('okay', 'The output directory seems correctly specified.')

    # Now we try to load the graph_table into a networkx.Graph() with weighted edges.
    """
    from: https://networkx.github.io/documentation/stable/reference/classes/index.html#which-graph-class-should-i-use
    Networkx Class      Type            Self-loops allowed      Parallel edges allowed
    Graph               undirected      Yes                     No
    DiGraph             directed        Yes                     No
    MultiGraph          undirected      Yes                     Yes
    MultiDiGraph        directed        Yes                     Yes
    """
    start_time = time.time()
    graph = net.read_edgelist(graph_table, nodetype=str, data=(('weight',float),), create_using=net.Graph())
    consoleout('okay', f'{parameters.graph} converted into a Graph object in {time.time() - start_time} s.')



    # 1st PART:
    # Here we want to fill the vConTACT2 csv output with the metas (taxonomy, etc) provided by INPHARED:
    # load csv pandas dataframe:
    start_time = time.time()
    csv = pnd.read_csv(csv_table, header = 0)
    consoleout('okay', f'{parameters.csv} loaded into memory in {time.time() - start_time} s.')

    # load metas pandas dataframe:
    start_time = time.time()
    metas = pnd.read_csv(metas_table, header = 0, sep='\t')
    consoleout('okay', f'{parameters.metas} loaded into memory in {time.time() - start_time} s.')

    # get the list of vOTUs:
    start_time = time.time()
    # These should be all our vOTUs (just a list of their names). These are the vOTUs contained in the "csv".
    # This should be equivalent to 'grep -Eo "vOTU_.*?," genome_by_genome_overview.csv | sort | uniq | wc -l'
    votus = csv[csv['Genome'].str.startswith(parameters.prefix)]['Genome'].tolist()
    # These should be all vOTUs contained in the graph:
    # This should be equivalent to 'grep -Eo "vOTU_.*? " c1.ntw | sort | uniq | wc -l'
    votus_ingraph = [votu for votu in votus if graph.has_node(votu) == True ]
    consoleout("okay", f"List of {len(votus_ingraph)}/{len(votus)} vOTU contained in the graph got in {time.time() - start_time} s.")

    # fill the vConTACT2 csv output with the taxonomy provided by INPHARED:
    start_time = time.time()
    csv_edit = fillWithMetas(csv, metas)
    consoleout('okay', f'{parameters.csv} updated with INPHARED db in {time.time() - start_time} s.')
    csv_edit.to_excel(parameters.output + 'csv_edit_' + parameters.suffix + '.xlsx' )# save the table!


    
    # 2nd PART:
    # Run the main algorithm. For every viral scaffold, get the most probable taxonomy: 
    start_time = time.time()
    df_results = clusterExtractor(graph, csv_edit, parameters.output, parameters.suffix, parameters.prefix)
    nC = len(df_results[df_results["Level"].str.startswith('C')])
    nN = len(df_results[df_results["Level"].str.startswith('N')])
    nG = len(df_results[df_results["Level"].str.startswith('G')])
    nF = len(df_results[df_results["Level"].str.startswith('F')])
    nA = len(df_results[df_results["Level"].str.startswith('A')])
    consoleout("okay", f'Taxonomy of {len(votus)} vOTUs (C:{nC}; N:{nN}; G:{nG}; F:{nF}; A:{nA}) obtained in {time.time() - start_time} s.')
    

    
    # 3rd PART:
    # generate the intractive subgraphs:
    # Edgse' color and width are relative to the weight. So here we first need to compute the max weight.
    # This way, colors and widths will be RELATIVE to the current maximum. In order to make subgraphs
    # comparable between different graphanalyzer.py calls, the max_weight need to be constrained.
    # After several tests, it appears to be already constrained to 300 by vConTACT2. 
    graph_table.close() # It was still opened for the networkx.Graph() creation.
    graph_table = open(parameters.graph, 'r')
    arrows = pnd.read_csv(graph_table, header = None, sep=' ')
    global_weights = list(arrows[2]) # get all weights in a list
    max_weight = 300
    if max(global_weights) > max_weight:
        consoleout("warning", f"Max weight here is > {max_weight} ({max(global_weights)}). {max_weight} will be used anyway. Please contact the developer.")
    consoleout("okay", f"Starting to generate {len(votus_ingraph)} subraphs with {parameters.threads} threads.")
    start_time = time.time()
    subgraphCreator(graph, csv_edit, df_results, parameters.output, parameters.suffix, max_weight, parameters.prefix, parameters.threads) 
    consoleout("okay", f'Generated {len(votus_ingraph)} interactive subgraphs in {time.time() - start_time} s.')
    

    print(f"Ending graphanalyzer.py v{version} on {datetime.datetime.now()}")
