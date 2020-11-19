# Coded by Gioele Lazzari (gioele.lazza@studenti.univr.it)
software = "db_manager.py"
version = "0.1.1"

import sys, os, argparse, logging, subprocess



parser = argparse.ArgumentParser(
    description = 'This script is designed to manage all the databases used '
                  'in the MetaPhage nextflow pipeline. This script is intended to be '
                  'placed in MetaPhage/bin.',
    formatter_class = argparse.RawTextHelpFormatter)

options = parser.add_argument_group("Options")

options.add_argument('-v', '--version', action='version', version= software + " v" + version)
# what follow are ALL the parameters passed by nextflow
options.add_argument('-p', '--mod_phix', dest='mod_phix', metavar='STRING', default=None)
options.add_argument('-p1', '--file_phix_alone', dest='file_phix_alone', metavar='PATH', default=None)
options.add_argument('-k', '--mod_kraken2', dest='mod_kraken2', metavar='STRING', default=None)
options.add_argument('-k1', '--file_kraken2_db', dest='file_kraken2_db', metavar='PATH', default=None)
options.add_argument('-b', '--mod_vibrant', dest='mod_vibrant', metavar='STRING', default=None)
options.add_argument('-b1', '--file_vibrant_db', dest='file_vibrant_db', metavar='PATH', default=None)
options.add_argument('-g', '--mod_phigaro', dest='mod_phigaro', metavar='STRING', default=None)
options.add_argument('-g1', '--file_phigaro_config', dest='file_phigaro_config', metavar='PATH', default=None)
options.add_argument('-s', '--mod_virsorter2', dest='mod_virsorter2', metavar='STRING', default=None)
options.add_argument('-s1', '--file_virsorter2_db', dest='file_virsorter2_db', metavar='PATH', default=None)



def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)



def manage(projectDir, 
           mod_phix, file_phix_alone,
           mod_kraken2, file_kraken2_db,
           mod_vibrant, file_vibrant_db,
           mod_phigaro, file_phigaro_config,
           mod_virsorter2, file_virsorter2_db):

    
    def checkCreate(path):
        if not os.path.exists(path):
            os.makedirs(path)
    
    def makeChannel(asset, content):
        # this file is the one actually opened in the process
        f = open(varDir + asset, "w")
        f.write(content)
        f.close()
        # this file in the working dir is just to temporize the process
        f = open(asset, "w")
        f.write(content)
        f.close()

    def echo(phrase):
        if phrase == "\n":
            sys.stdout.write("\n\n")
        else:
            sys.stdout.write("db_manager: " + phrase + "\n")


    
    varDir = projectDir + "bin/groovy_vars/"
    checkCreate(varDir)
    echo("\n")


    #################################
    # phix ##########################
    #################################
    rel_path = "db/phix/"
    mod_folder = projectDir + rel_path
    checkCreate(mod_folder)
    if mod_phix == "custom":
        if file_phix_alone == "-":
            error('With --mod_phix custom you have to specify also --file_phix_alone')
        else:
            makeChannel("file_phix_alone", file_phix_alone)

    elif mod_phix == "phiX174":
        if not os.path.exists(mod_folder + "phiX174.fasta"):
            url = "https://www.ebi.ac.uk/ena/browser/api/fasta/AF176027.1?download=true"
            echo("Downloading " + mod_folder + "phiX174.fasta" + " ...")
            os.popen('wget -O ' + mod_folder + "phiX174.fasta" + ' ' + url).read()
            echo("OK")
        else:
            echo(mod_folder + "phiX174.fasta" + ' already present!')
        makeChannel("file_phix_alone", "db/phix/phiX174.fasta")

    elif mod_phix == "WA11":
        if not os.path.exists(mod_folder + "WA11.fasta"):
            url = "https://www.ebi.ac.uk/ena/browser/api/fasta/DQ079895.1?download=true"
            echo("Downloading " + mod_folder + "WA11.fasta" + " ...")
            os.popen('wget -O ' + mod_folder + "WA11.fasta" + ' ' + url).read() 
            echo("OK")
        else:
            echo(mod_folder + "WA11.fasta" + ' already present!')
        makeChannel("file_phix_alone", "db/phix/WA11.fasta")

    
    #################################
    # kraken2 #######################
    #################################
    if mod_kraken2 == "custom":
        if file_kraken2_db == "-":
            error('With --mod_kraken2 custom you have to specify also --file_kraken2_db')
        else:
            makeChannel("file_kraken2_db", file_kraken2_db)
    
    elif mod_kraken2 == "miniBAV":
        rel_path = "db/kraken2/miniBAV/"
        mod_folder = projectDir + rel_path
        checkCreate(mod_folder)
        if not os.path.exists(mod_folder + "hash.k2d") or not os.path.exists(mod_folder + "opts.k2d") or not os.path.exists(mod_folder + "taxo.k2d"):
            url = "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz"
            echo("Downloading " + mod_folder + "*" + " ...")
            os.popen('wget -O ' + mod_folder + "minikraken2_v1_8GB_201904.tgz" + ' ' + url).read() # note for future: without read(), wget won't pass  to tar!!!
            echo("Decompressing ...")
            os.popen('tar zxvf %s --directory %s' % (mod_folder + "minikraken2_v1_8GB_201904.tgz", mod_folder)).read() # note for future: without read(), tar won't extract anything!!!
            os.popen('mv %s %s' % (mod_folder + "minikraken2_v1_8GB/*", mod_folder)).read()
            os.popen('rm -rf %s' % (mod_folder + "minikraken2_v1_8GB")).read()
            os.popen('rm %s' % (mod_folder + "minikraken2_v1_8GB_201904.tgz")).read()
            echo("OK") 
        else:
            echo(mod_folder + "*" + ' already present!')
        makeChannel("file_kraken2_db", rel_path)
    
    elif mod_kraken2 == "miniBAVH":
        rel_path = "db/kraken2/miniBAVH/"
        mod_folder = projectDir + rel_path
        checkCreate(mod_folder)
        if not os.path.exists(mod_folder + "hash.k2d") or not os.path.exists(mod_folder + "opts.k2d") or not os.path.exists(mod_folder + "taxo.k2d"):
            url = "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz"
            echo("Downloading " + mod_folder + "*" + " ...")
            os.popen('wget -O ' + mod_folder + "minikraken2_v2_8GB_201904.tgz" + ' ' + url).read() # note for future: without read(), wget won't pass  to tar!!!
            echo("Decompressing ...")
            os.popen('tar zxvf %s --directory %s' % (mod_folder + "minikraken2_v2_8GB_201904.tgz", mod_folder)).read() # note for future: without output.read(), tar won't extract anything!!!
            os.popen('mv %s %s' % (mod_folder + "minikraken2_v2_8GB_201904_UPDATE/*", mod_folder)).read()
            os.popen('rm -rf %s' % (mod_folder + "minikraken2_v2_8GB_201904_UPDATE")).read()
            os.popen('rm %s' % (mod_folder + "minikraken2_v2_8GB_201904.tgz")).read()
            echo("OK")
        else:
            echo(mod_folder + "*" + ' already present!')
        makeChannel("file_kraken2_db", rel_path)



    #################################
    # vibrant #######################
    #################################
    if mod_vibrant == "custom":
        if file_vibrant_db == "-":
            error('With --mod_vibrant custom you have to specify also --file_vibrant_db')
        else:
            makeChannel("file_vibrant_db", file_vibrant_db)
    
    elif mod_vibrant == "standard":
        rel_path = "db/vibrant/standard/"
        mod_folder = projectDir + rel_path
        checkCreate(mod_folder)
        checkCreate(mod_folder + "databases")
        checkCreate(mod_folder + "files")
        if not os.path.exists(mod_folder + "databases/KEGG_profiles_prokaryotes.HMM.h3p") or not os.path.exists(mod_folder + "databases/Pfam-A_v32.HMM.h3p") or not os.path.exists(mod_folder + "databases/VOGDB94_phage.HMM.h3p"):
            echo("Downloading " + mod_folder + "*" + " ...")
            os.popen('wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz').read()
            os.popen('wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz').read()
            os.popen('wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz').read()
            echo("Decompressing ...")
            os.popen('mkdir profile_names').read()
            # vog
            os.popen('tar -xzf vog.hmm.tar.gz').read()
            os.popen('for v in VOG*.hmm; do cat $v >> vog_temp.HMM; done').read()
            os.popen('rm VOG0*.hmm').read()
            os.popen('rm VOG1*.hmm').read()
            os.popen('rm VOG2*.hmm').read()
            os.popen('wget -O profile_names/VIBRANT_vog_profiles.txt https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/databases/profile_names/VIBRANT_vog_profiles.txt').read()
            os.popen('hmmfetch -o VOGDB94_phage.HMM -f vog_temp.HMM profile_names/VIBRANT_vog_profiles.txt').read()
            os.popen('hmmpress VOGDB94_phage.HMM').read()
            os.popen('rm vog_temp.HMM vog.hmm.tar.gz').read()
            # pfam
            os.popen('gunzip Pfam-A.hmm.gz').read()
            os.popen('mv Pfam-A.hmm Pfam-A_v32.HMM').read()
            os.popen('hmmpress Pfam-A_v32.HMM').read()
            # kegg
            os.popen('tar -xzf profiles.tar.gz').read()
            os.popen('for k in profiles/K*.hmm; do cat $k >> kegg_temp.HMM; done').read()
            os.popen('wget -O profile_names/VIBRANT_kegg_profiles.txt https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/databases/profile_names/VIBRANT_kegg_profiles.txt').read()
            os.popen('hmmfetch -o KEGG_profiles_prokaryotes.HMM -f kegg_temp.HMM profile_names/VIBRANT_kegg_profiles.txt').read()
            os.popen('hmmpress KEGG_profiles_prokaryotes.HMM').read()
            os.popen('rm -R profiles').read()
            os.popen('rm kegg_temp.HMM profiles.tar.gz').read()
            # move to mod_folder/databases/
            os.popen('rm -R profile_names').read()
            os.popen('mv VOG* %sdatabases/' % (mod_folder)).read()
            os.popen('mv Pfam* %sdatabases/' % (mod_folder)).read()
            os.popen('mv KEGG* %sdatabases/' % (mod_folder)).read()
            # now fill mod_folder/files/
            os.popen('wget -O %sfiles/VIBRANT_AMGs.tsv https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/files/VIBRANT_AMGs.tsv' % (mod_folder)).read()
            os.popen('wget -O %sfiles/VIBRANT_KEGG_pathways_summary.tsv https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/files/VIBRANT_KEGG_pathways_summary.tsv' % (mod_folder)).read()
            os.popen('wget -O %sfiles/VIBRANT_categories.tsv https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/files/VIBRANT_categories.tsv' % (mod_folder)).read()
            os.popen('wget -O %sfiles/VIBRANT_machine_model.sav https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/files/VIBRANT_machine_model.sav' % (mod_folder)).read()
            os.popen('wget -O %sfiles/VIBRANT_names.tsv https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/files/VIBRANT_names.tsv' % (mod_folder)).read()
            
            echo("OK") 
        else:
            echo(mod_folder + "*" + ' already present!')
        makeChannel("file_vibrant_db", rel_path)

    elif mod_vibrant == "legacy":
        rel_path = "db/vibrant/legacy/"
        mod_folder = projectDir + rel_path
        checkCreate(mod_folder)
        if not os.path.exists(mod_folder + "KEGG_profiles_prokaryotes.HMM") or not os.path.exists(mod_folder + "Pfam-A_v32.HMM") or not os.path.exists(mod_folder + "VOGDB94_phage.HMM"):
            echo("Downloading " + mod_folder + "*" + " ...")
            os.popen('wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz').read()
            os.popen('wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz').read()
            os.popen('wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz').read()
            echo("Decompressing ...")
            os.popen('mkdir profile_names').read()
            # vog
            os.popen('tar -xzf vog.hmm.tar.gz').read()
            os.popen('for v in VOG*.hmm; do cat $v >> vog_temp.HMM; done').read()
            os.popen('rm VOG0*.hmm').read()
            os.popen('rm VOG1*.hmm').read()
            os.popen('rm VOG2*.hmm').read()
            os.popen('wget -O profile_names/VIBRANT_vog_profiles.txt https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/databases/profile_names/VIBRANT_vog_profiles.txt').read()
            os.popen('hmmfetch -o VOGDB94_phage.HMM -f vog_temp.HMM profile_names/VIBRANT_vog_profiles.txt').read()
            os.popen('hmmpress VOGDB94_phage.HMM').read()
            os.popen('rm vog_temp.HMM vog.hmm.tar.gz').read()
            # pfam
            os.popen('gunzip Pfam-A.hmm.gz').read()
            os.popen('wget -O profile_names/VIBRANT_pfam-plasmid_profiles.txt https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/7a59d65ee21e7a23f9579914b9eb3bc83149b62f/databases/profile_names/VIBRANT_pfam-plasmid_profiles.txt').read()
            os.popen('hmmfetch -o Pfam-A_plasmid_v32.HMM -f Pfam-A.hmm profile_names/VIBRANT_pfam-plasmid_profiles.txt').read()
            os.popen('wget -O profile_names/VIBRANT_pfam-phage_profiles.txt https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/7a59d65ee21e7a23f9579914b9eb3bc83149b62f/databases/profile_names/VIBRANT_pfam-phage_profiles.txt').read()
            os.popen('hmmfetch -o Pfam-A_phage_v32.HMM -f Pfam-A.hmm profile_names/VIBRANT_pfam-phage_profiles.txt').read()
            os.popen('mv Pfam-A.hmm Pfam-A_v32.HMM').read()
            os.popen('hmmpress Pfam-A_plasmid_v32.HMM').read()
            os.popen('hmmpress Pfam-A_phage_v32.HMM').read()
            os.popen('hmmpress Pfam-A_v32.HMM').read()
            # kegg
            os.popen('tar -xzf profiles.tar.gz').read()
            os.popen('for k in profiles/K*.hmm; do cat $k >> kegg_temp.HMM; done').read()
            os.popen('wget -O profile_names/VIBRANT_kegg_profiles.txt https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/master/databases/profile_names/VIBRANT_kegg_profiles.txt').read()
            os.popen('hmmfetch -o KEGG_profiles_prokaryotes.HMM -f kegg_temp.HMM profile_names/VIBRANT_kegg_profiles.txt').read()
            os.popen('hmmpress KEGG_profiles_prokaryotes.HMM').read()
            os.popen('rm -R profiles').read()
            os.popen('rm kegg_temp.HMM profiles.tar.gz').read()
            # move to mod_folder/databases/
            os.popen('rm -R profile_names').read()
            os.popen('mv VOG* %s' % (mod_folder)).read()
            os.popen('mv Pfam* %s' % (mod_folder)).read()
            os.popen('mv KEGG* %s' % (mod_folder)).read()
            # now fill mod_folder/files/
            os.popen('wget -O %sVIBRANT_AMGs.tsv https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/7a59d65ee21e7a23f9579914b9eb3bc83149b62f/files/VIBRANT_AMGs.tsv' % (mod_folder)).read()
            os.popen('wget -O %sVIBRANT_KEGG_pathways_summary.tsv https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/7a59d65ee21e7a23f9579914b9eb3bc83149b62f/files/VIBRANT_KEGG_pathways_summary.tsv' % (mod_folder)).read()
            os.popen('wget -O %sVIBRANT_categories.tsv https://raw.githubusercontent.com/AnantharamanLab/VIBRANT/7a59d65ee21e7a23f9579914b9eb3bc83149b62f/files/VIBRANT_categories.tsv' % (mod_folder)).read()
            os.popen('wget -O %sVIBRANT_machine_model.sav https://github.com/AnantharamanLab/VIBRANT/raw/7a59d65ee21e7a23f9579914b9eb3bc83149b62f/files/VIBRANT_machine_model.sav' % (mod_folder)).read()
            os.popen('wget -O %sVIBRANT_names.tsv https://github.com/AnantharamanLab/VIBRANT/raw/7a59d65ee21e7a23f9579914b9eb3bc83149b62f/files/VIBRANT_names.tsv' % (mod_folder)).read()
            
            echo("OK") 
        else:
            echo(mod_folder + "*" + ' already present!')
        makeChannel("file_vibrant_db", rel_path)


    #################################
    # phigaro #######################
    #################################
    if mod_phigaro == "custom":
        if file_phigaro_config == "-":
            error('With --mod_phigaro custom you have to specify also --file_phigaro_config')
        else:
            makeChannel("file_phigaro_config", file_phigaro_config)
    
    elif mod_phigaro == "standard":
        rel_path = "db/phigaro/standard/"
        mod_folder = projectDir + rel_path
        checkCreate(mod_folder)
        if not os.path.exists(mod_folder + "allpvoghmms.h3f") or not os.path.exists(mod_folder + "allpvoghmms.h3i") or not os.path.exists(mod_folder + "allpvoghmms.h3m") or not os.path.exists(mod_folder + "allpvoghmms.h3p"):
            echo("Downloading " + mod_folder + "*" + " ...")
            stri = os.popen('wget -O %sallpvoghmms http://download.ripcm.com/phigaro/allpvoghmms' % (mod_folder)).read()
            os.popen('wget -O %sallpvoghmms.h3f http://download.ripcm.com/phigaro/allpvoghmms.h3f' % (mod_folder)).read()
            os.popen('wget -O %sallpvoghmms.h3i http://download.ripcm.com/phigaro/allpvoghmms.h3i' % (mod_folder)).read()
            os.popen('wget -O %sallpvoghmms.h3m http://download.ripcm.com/phigaro/allpvoghmms.h3m' % (mod_folder)).read()
            os.popen('wget -O %sallpvoghmms.h3p http://download.ripcm.com/phigaro/allpvoghmms.h3p' % (mod_folder)).read()
            
            echo("OK") 
        else:
            echo(mod_folder + "*" + ' already present!')
        makeChannel("file_phigaro_config", rel_path)


    
    #################################
    # virsorter2 ####################
    #################################
    if mod_virsorter2 == "custom":
        if file_virsorter2_db == "-":
            error('With --mod_virsorter2 custom you have to specify also --file_virsorter2_db')
        else:
            makeChannel("file_virsorter2_db", file_virsorter2_db)
    
    elif mod_virsorter2 == "standard":
        rel_path = "db/virsorter2/standard/"
        mod_folder = projectDir + rel_path
        checkCreate(mod_folder)
        if not os.path.exists(mod_folder + "file1") or not os.path.exists(mod_folder + "file2"):
            echo("Downloading " + mod_folder + "*" + " ...")
            
            echo("OK") 
        else:
            echo(mod_folder + "*" + ' already present!')
        makeChannel("file_virsorter2_db", rel_path)
    



if __name__ == "__main__":

    parameters = parser.parse_args()

    # understand project directory
    projectDir = os.path.realpath(__file__).replace("bin/db_manager.py", "")

    # core function
    manage(projectDir, 
           parameters.mod_phix, parameters.file_phix_alone,
           parameters.mod_kraken2, parameters.file_kraken2_db,
           parameters.mod_vibrant, parameters.file_vibrant_db,
           parameters.mod_phigaro, parameters.file_phigaro_config,
           parameters.mod_virsorter2, parameters.file_virsorter2_db)

    
