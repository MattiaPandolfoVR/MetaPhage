#!/usr/bin/env python3
# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
software = "phigaro_config_creator.py"
version = "0.0.1"

import  os, sys

db_path = sys.argv[1]

if __name__ == "__main__":

    path_prodigal = os.popen('which prodigal').read()
    path_hmmsearch = os.popen('which hmmsearch').read()
    
    #projectDir = os.path.realpath(__file__).replace("bin/python/phigaro_config_creator.py", "")
    #path_allpvoghmms = projectDir + "db/phigaro/standard/allpvoghmms"
    if os.path.exists(db_path):
      database_dir = os.path.realpath(db_path)
      path_allpvoghmms = database_dir + "/phigaro/standard/allpvoghmms"
    else:
      exit()

    file_config = open("config.yml", "w")

    file_config.write("""
hmmer:
  bin: %s
  e_value_threshold: 0.00445
  pvog_path: %s
phigaro:
  mean_gc: 0.46354823199323625
  penalty_black: 2.2
  penalty_white: 0.7
  threshold_max_abs: 52.96
  threshold_max_basic: 46.0
  threshold_max_without_gc: 11.42
  threshold_min_abs: 50.32
  threshold_min_basic: 45.39
  threshold_min_without_gc: 11.28
  window_len: 32
prodigal:
  bin: %s
""" % (path_hmmsearch, path_allpvoghmms, path_prodigal))
    
    file_config.close()
