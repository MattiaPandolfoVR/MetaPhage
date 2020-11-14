# Coded by Gioele Lazzari (gioele.lazza@studenti.univr.it)
software = "db_manager.py"
version = "0.0.3"

import sys, os, argparse, logging, subprocess


#############
# ARGPARSER #
#############

parser = argparse.ArgumentParser(
    description = 'This script is designed to manage all the databases used '
                  'in the MetaPhage nextflow pipeline. This script is intended to be '
                  'placed in MetaPhage/bin. This script works in conjunction with '
                  '2 files, "file_table.csv" and "last_settings.csv", both placed in MetaPhage/db.',
    formatter_class = argparse.RawTextHelpFormatter)

options = parser.add_argument_group("Options")

options.add_argument('-v', '--version', action='version', version= software + " v" + version)
# what follow are ALL the parameters passed by nextflow
# phix
options.add_argument('-p', '--mod_phix', dest='mod_phix', metavar='STRING', default=None)
options.add_argument('-p1', '--file_phix_alone', dest='file_phix_alone', metavar='PATH', default=None)
# kraken2
options.add_argument('-k', '--mod_kraken2', dest='mod_kraken2', metavar='STRING', default=None)
options.add_argument('-k1', '--file_kraken2_db', dest='file_kraken2_db', metavar='PATH', default=None)


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


##################
# CORE FUNCTIONS #
##################


def manage(projectDir, df_file_table, df_last_settings,
           # phix
           mod_phix, file_phix_alone,
           # kraken2
           mod_kraken2, file_kraken2_db):

    # PART 1: update last preferences

    #################################
    # phix ##########################
    #################################
    if (mod_phix != "-"):
        # subset dataframe for the given program
        rows_phix = df_file_table[df_file_table['file'].str.contains('file_phix')] 
        # get every modality for that program
        mod_list = list(set(rows_phix['modality'].tolist()))
        mod_list.append("custom") 
        # get every requested file for that program
        file_list = list(set(rows_phix['file'].tolist())) 

        # if the user specified modality is allowed, save this run's parameters into df_last_settings
        if mod_phix in mod_list: 
            if mod_phix == "custom" and file_phix_alone == "-":
                error('For a program with its databses/models in "custom" mode you have to specify all thier file paths')
            else:
                for item in file_list:
                    df_last_settings.loc[df_last_settings["file"] == item, ["modality"]] = mod_phix

                    if mod_phix != "custom":
                        path_to_use = df_file_table.loc[(df_file_table['file']== item) & (df_file_table['modality'] == mod_phix), 'filepath'].values[0]
                        df_last_settings.loc[df_last_settings["file"] == item, ["filepath"]] = path_to_use

                    elif mod_phix == "custom" and file_phix_alone != "-":
                        df_last_settings.loc[df_last_settings["file"] == item, ["filepath"]] = file_phix_alone                    
        else:
            error("Modality not allowed. Allowed modalities are " + str(mod_list))

    #################################
    # kraken2 #######################
    #################################
    if (mod_kraken2 != "-"):
        # subset dataframe for the given program
        rows_kraken2 = df_file_table[df_file_table['file'].str.contains('file_kraken2')] 
        # get every modality for that program
        mod_list = list(set(rows_kraken2['modality'].tolist()))
        mod_list.append("custom") 
        # get every requested file for that program
        file_list = list(set(rows_kraken2['file'].tolist())) 

        # if the user specified modality is allowed, save this run's parameters into df_last_settings
        if mod_kraken2 in mod_list: 
            if mod_kraken2 == "custom" and file_kraken2_db == "-":
                error('For a program with its databses/models in "custom" mode you have to specify all thier file paths')
            else:
                for item in file_list:
                    df_last_settings.loc[df_last_settings["file"] == item, ["modality"]] = mod_kraken2

                    if mod_kraken2 != "custom":
                        path_to_use = df_file_table.loc[(df_file_table['file']== item) & (df_file_table['modality'] == mod_kraken2), 'filepath'].values[0]
                        df_last_settings.loc[df_last_settings["file"] == item, ["filepath"]] = path_to_use

                    elif mod_kraken2 == "custom" and file_kraken2_db != "-":
                        df_last_settings.loc[df_last_settings["file"] == item, ["filepath"]] = file_kraken2_db                    
        else:
            error("Modality not allowed. Allowed modalities are " + str(mod_list))


    # PART 2: save current preferences
    df_last_settings.to_csv(projectDir + "db/last_settings.csv", index = False, sep = ';')


    # PART 3: export current preferences in individual file readable by groovy
    # create a subfolder for containg all the variables
    if not os.path.exists(projectDir + "db/groovy_vars"):
        os.makedirs(projectDir + "db/groovy_vars")
    for row in df_last_settings.itertuples():
        # this file is the one actually opened in the process
        f= open(projectDir + "db/groovy_vars/" + row.file, "w")
        f.write(row.filepath)
        f.close()
        # this file in the working dir is just to temporize the process
        f= open(row.file, "w")
        f.write(row.filepath)
        f.close()


    # PART 4: check file presence and download it if absent
    logging.info("Current settings are:\n\n" + df_last_settings.to_string() + "\n\n\n\n\n")
    
    #################################
    # phix ##########################
    #################################
    # subset dataframe for the given program
    rows_phix = df_last_settings[df_last_settings['file'].str.contains('file_phix')]
    # for every file required by the program:
    for row in rows_phix.itertuples():

        if row.modality == "custom":
            if os.path.exists(projectDir + row.filepath) == False:
                error("Can't find the file you specified: " + projectDir + row.filepath )
        # if the program is NOT set in "custom" mode, manage its files automatically:
        else:
            # extract the predefined path of this file
            relative_path = df_file_table.loc[(df_file_table['file'].str.contains('file_phix')) & 
                                              (df_file_table['modality'] == row.modality), 'filepath'].values[0]
            path_to_check = projectDir + relative_path
            # if this file is NOT already downloaded:
            if os.path.exists(path_to_check) == False:
                # extract the predefined url for this file to be downloaded
                url = df_file_table.loc[(df_file_table['file'].str.contains('file_phix')) & 
                                        (df_file_table['modality'] == row.modality), 'url'].values[0]
                # download the file with wget
                # as described here https://janakiev.com/blog/python-shell-commands/ 
                logging.info('%s is missing. It will downloaded from %s with the following command:\n\nwget -O %s %s \n'
                             %(relative_path, url, path_to_check, url))
                # launch the command and get its verbose into log
                process = subprocess.Popen(['wget', '-O', path_to_check, url], 
                           stdout=subprocess.PIPE,  stderr=subprocess.STDOUT)
                def check_io():
                    while True:
                        output = process.stdout.readline().decode()
                        if output:
                            logging.info(output.rstrip("\n")) # rstrip remove the trailing endline
                        else:
                            break
                # keep checking stdout/stderr until the child exits
                while process.poll() is None:
                    check_io()
                logging.info("This download is complete.\n\n\n\n\n")
                
            else: # this file IS already downloaded in the specified location
                logging.info('%s is already downloaded in the specified path. No need to download it again.\n\n\n\n\n'
                             %(relative_path))

    #################################
    # kraken2 #######################
    #################################
    # subset dataframe for the given program
    rows_kraken2 = df_last_settings[df_last_settings['file'].str.contains('file_kraken2')]
    # for every file required by the program:
    for row in rows_kraken2.itertuples():

        if row.modality == "custom":
            if os.path.exists(projectDir + row.filepath) == False:
                error("Can't find the file you specified: " + projectDir + row.filepath )
        # if the program is NOT set in "custom" mode, manage its files automatically:
        else:
            # extract the predefined path of this file
            relative_path = df_file_table.loc[(df_file_table['file'].str.contains('file_kraken2')) & 
                                              (df_file_table['modality'] == row.modality), 'filepath'].values[0]
            path_to_check = projectDir + relative_path
            # if this file is NOT already downloaded:
            if os.path.exists(path_to_check) == False:
                # extract the predefined url for this file to be downloaded
                url = df_file_table.loc[(df_file_table['file'].str.contains('file_kraken2')) & 
                                        (df_file_table['modality'] == row.modality), 'url'].values[0]
                # download the file with wget
                # as described here https://janakiev.com/blog/python-shell-commands/ 
                logging.info('%s is missing. It will downloaded from %s with the following command:\n\nwget -O %s %s \n'
                             %(relative_path, url, path_to_check.replace("hash.k2d", "archive.tgz"), url))
                # launch the command and get its verbose into log
                process = subprocess.Popen(['wget', '-O', path_to_check.replace("hash.k2d", "archive.tgz"), url], 
                           stdout=subprocess.PIPE,  stderr=subprocess.STDOUT)
                def check_io():
                    while True:
                        output = process.stdout.readline().decode()
                        if output:
                            logging.info(output.rstrip("\n")) # rstrip remove the trailing endline
                        else:
                            break
                # keep checking stdout/stderr until the child exits
                while process.poll() is None:
                    check_io()
                logging.info("This download is complete.\n\n\n\n\n")

                # extract pre-compiled databases
                logging.info("Decompressing %s...\n\n\n\n\n" % (path_to_check.replace("hash.k2d", "archive.tgz")))
                os.popen('tar zxvf %s' % (path_to_check.replace("hash.k2d", "archive.tgz")))
                os.popen('mv %s %s' % (path_to_check.replace("hash.k2d", "minikraken2_v*_8GB/*"), path_to_check.replace("hash.k2d", "")))
                os.popen('rm -r %s' % (path_to_check.replace("hash.k2d", "minikraken2_v*_8GB")))
                os.popen('rm %s' % (path_to_check.replace("hash.k2d", "archive.tgz")))
                
            else: # this file IS already downloaded in the specified location
                logging.info('%s is already downloaded in the specified path. No need to download it again.\n\n\n\n\n'
                             %(relative_path))
 
if __name__ == "__main__":

    parameters = parser.parse_args()

    # create the .log file
    logging.basicConfig(filename ="./db_manager.log", filemode='w', level = logging.INFO, # level sets the threshold
                        format = '%(asctime)s %(levelname)s: %(message)s',
                        datefmt = '%H:%M:%S') 
    logging.info('Processing of databases is started.\n')

    # understand project dirrectory
    stream = os.popen('pwd')
    logging.info("workingDir: " + stream.read())
    projectDir = os.path.realpath(__file__).replace("bin/db_manager.py", "")
    logging.info("projectDir: " + projectDir + "\n")
    
    # convert system files to pandas dataframes
    df_file_table = pnd.read_csv(projectDir + "db/file_table.csv", header = 0, sep = ';')
    df_last_settings = pnd.read_csv(projectDir + "db/last_settings.csv", header = 0, sep = ';')

    # core function
    manage(projectDir, df_file_table, df_last_settings,
           # phix
           parameters.mod_phix, parameters.file_phix_alone,
           # kraken2
           parameters.mod_kraken2, parameters.file_kraken2_db)

    
