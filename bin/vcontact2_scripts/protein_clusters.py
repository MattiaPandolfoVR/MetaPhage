"""Protein_clusters.py"""

import os
import gzip
from Bio import SeqIO
import pandas as pd
import logging
import subprocess
import numpy as np


logger = logging.getLogger(__name__)


def merge_aa(user_aa_fp, ref_db_fp, merged_aa_fp):
    """

    :param user_aa_fp:
    :param ref_db_fp:
    :return:
    """

    with open(merged_aa_fp, 'w') as merged_aa_fh:
        for aa_fp in [user_aa_fp, ref_db_fp]:
            if '.gz' in aa_fp:
                with gzip.open(aa_fp, 'rt') as aa_fh:
                    seq_records = SeqIO.parse(aa_fh, 'fasta')
                    SeqIO.write(seq_records, merged_aa_fh, 'fasta')
            else:
                with open(aa_fp, 'rU') as aa_fh:
                    seq_records = SeqIO.parse(aa_fh, 'fasta')
                    SeqIO.write(seq_records, merged_aa_fh, 'fasta')

    return merged_aa_fp


def make_blast_db(aa_fp):

    """
    :param aa_fp: Amino acid fasta file file path
    :return: Output database file path
    """

    out_db = '{}.db'.format(aa_fp.rsplit('.', 1)[0])
    makeblastcmd = 'makeblastdb -in {} -input_type fasta -dbtype prot -hash_index -out {}'.format(aa_fp, out_db)

    logger.debug("Creating BLAST DB...")
    subprocess.check_call(makeblastcmd, shell=True)

    return out_db


def run_blastp(aa_fp, db_fp, evalue: float, cpu: int, blastp_out_fp):
    """
    :param aa_fp: File path for amino acid fast file
    :param db_fp: BLAST database file path
    :param evalue:
    :param cpu:
    :param blastp_out_fp:
    :return:
    """

    blastp_cmd = 'blastp -task blastp -query {} -db {} -out {} -evalue {} -outfmt 6 -num_threads {}'.format(
        aa_fp, db_fp, blastp_out_fp, evalue, cpu)

    logger.debug("Running BLASTP...")
    subprocess.check_call(blastp_cmd, shell=True)

    return blastp_out_fp


def make_diamond_db(aa_fp, db_dir, cpu: int):

    diamond_db_bp = os.path.join(db_dir, os.path.basename(aa_fp).rsplit('.', 1)[0])
    make_diamond_cmd = ['diamond', 'makedb', '--threads', str(cpu), '--in', aa_fp, '-d', diamond_db_bp]

    logger.info("Creating Diamond database...")
    res = subprocess.run(make_diamond_cmd, check=True, stdout=subprocess.PIPE)

    if res.returncode != 0:
        logger.error('Error creating Diamond database')
        exit(1)

    diamond_db_fp = diamond_db_bp + '.dmnd'

    return diamond_db_fp


def run_diamond(aa_fp, db_fp, cpu: int, diamond_out_fn):

    # More sensitive as an option?
    diamond_cmd = ['diamond', 'blastp', '--threads', str(cpu), '--sensitive', '-d', db_fp, '-q', aa_fp,
                   '-o', diamond_out_fn]

    logger.info("Running Diamond...")
    res = subprocess.run(diamond_cmd, check=True, stdout=subprocess.PIPE)

    if res.returncode != 0:
        logger.error('Error running Diamond')
        exit(1)

    return diamond_out_fn


def make_protein_clusters_mcl(blast_fp, out_p, inflation=2):
    """
    Args: 
        blast_fp (str): Path to blast results file
        inflation (float): MCL inflation value
        out_p (str): Output directory path
    Returns:
        str: fp for MCL clustering file
    """

    logger.debug("Generating abc file...")

    blast_fn = os.path.basename(blast_fp)
    abc_fn = '{}.abc'.format(blast_fn)
    abc_fp = os.path.join(out_p, abc_fn)
    subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(blast_fp, abc_fp), shell=True)

    logger.debug("Running MCL...")

    mci_fn = '{}.mci'.format(blast_fn)
    mci_fp = os.path.join(out_p, mci_fn)
    mcxload_fn = '{}_mcxload.tab'.format(blast_fn)
    mcxload_fp = os.path.join(out_p, mcxload_fn)
    subprocess.check_call("mcxload -abc {0} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {1}"
                          " -write-tab {2}".format(abc_fp, mci_fp, mcxload_fp), shell=True)

    mcl_clstr_fn = "{0}_mcl{1}.clusters".format(blast_fn, int(inflation*10))
    mcl_clstr_fp = os.path.join(out_p, mcl_clstr_fn)

    subprocess.check_call("mcl {0} -I {1} -use-tab {2} -o {3}".format(
        mci_fp, inflation, mcxload_fp, mcl_clstr_fp), shell=True)

    return mcl_clstr_fp


def make_protein_clusters_one(blast_fp, c1_bin, out_p, overlap: float, penalty: float, haircut: float):
    """
    Args:
        blast_fp (str): Blast results file
        out_p (str): output directory path
        overlap (int): hold
        penalty (int): hold
        haircut (int): hold
    Returns:
        str: fp for ClusterONE clustering file
    """

    # Grab first few columns of blastp output, effectively generating ClusterONE input
    blast_fn = os.path.basename(blast_fp)
    abc_fn = '{}.abc'.format(blast_fn)
    abc_fp = os.path.join(out_p, abc_fn)
    subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}.abc".format(blast_fp, abc_fp), shell=True)

    logger.debug("Running ClusterONE...")

    cluster_one_cmd = 'java -jar {} {}.abc --input-format edge_list --output-format csv ' \
                      '--max-overlap {} --penalty {} --haircut {}'.format(c1_bin, abc_fp, overlap, penalty, haircut)

    cluster_one_fn = "{}_one_{}_{}_{}.clusters".format(blast_fn, overlap, penalty, haircut)
    cluster_one_fp = os.path.join(out_p, cluster_one_fn)
    cluster_one_cmd += ' > {}'.format(cluster_one_fp)

    subprocess.check_call(cluster_one_cmd, shell=True)

    return cluster_one_fp


def build_clusters(fp, gene2genome, mode='ClusterONE'):
    """
        Build clusters given clusters file

        Args:
            fp (str): filepath of clusters file
            gene2genome (dataframe): A dataframe giving the protein and its genome.
            mode (str): clustering method
        Returns:
            tuple: dataframe of proteins, clusters, profiles and contigs
        """

    # Read MCL
    if mode == 'ClusterONE':
        clusters_df, name, c = load_one_clusters(fp)
    elif mode == 'MCL':
        clusters_df, name, c = load_mcl_clusters(fp)
    else:
        clusters_df, name, c = False
        logger.error("A mode must be selected. Use ClusterONE or MCL to generate PCs.")

    gene2genomeTMP=gene2genome.copy(deep=True)
    
    # Assign each prot to its cluster
    gene2genomeTMP.set_index("protein_id", inplace=True)  # id, contig, keywords, cluster

    #print("****** in module all'inizio")
    #print(gene2genomeTMP.head())
    #print("****************")

    
    x=0
    flat_list_prots = []
    flat_list_clust = []
    for prots in c:
        clust=name[x]
        x=x+1
        for item in prots:
            flat_list_prots.append(item)
            flat_list_clust.append(clust)
    data = {'protein_id':flat_list_prots,
            'cluster':flat_list_clust}
    df=pd.DataFrame(data)
    
    df.set_index("protein_id", inplace=True)
    
    left_merge=pd.merge(gene2genomeTMP,df,on="protein_id", how="left")
    #print(left_merge.to_string())
    gene2genomeTMP=left_merge

    #print("****** in module alla fine 2")
    #print(gene2genomeTMP.head())
    #print("****************")
    
    clust_num = []
    clust_list = []
    for clust, prots in gene2genomeTMP.groupby("cluster"):
         clust_num.append(prots.keywords.count())
         clust_list.append(clust)
         
    data = {'pc_id':clust_list,
            'annotated':clust_num}
    
    df_annotation=pd.DataFrame(data)
    df_annotation.set_index("pc_id", inplace=True)
    #gene2genome.reset_index(inplace=True)

    #print(clusters_df.head())
    #exit;
    
    left_merge=pd.merge(clusters_df,df_annotation,on="pc_id", how="left")
    clusters_df=left_merge.copy(deep=True)
    #clusters_df['annotated'] = clusters_df['annotated'].astype('int')
    #gene2genome.set_index("protein_id", inplace=True)
    
    for clust, prots in gene2genomeTMP.groupby("cluster"):
        if prots.keywords.count():
            keys = ";".join(prots.keywords.dropna().values).split(";")
            key_count = {}
            for k in keys:
                k = k.strip()
                try:
                    key_count[k] += 1
                except KeyError:
                    key_count[k] = 1

            clusters_df.loc[clust, "keys"] = "; ".join(["{} ({})".format(x, y) for x, y in key_count.items()])

    gene2genomeTMP.reset_index(inplace=True)
    clusters_df.reset_index(inplace=True)
    profiles_df = gene2genomeTMP.loc[:, ["contig_id", "cluster"]].drop_duplicates()
    profiles_df.columns = ["contig_id", "pc_id"]

    contigs_df = pd.DataFrame(gene2genomeTMP.fillna(0).groupby("contig_id").count().protein_id)
    contigs_df.index.name = "contig_id"
    contigs_df.columns = ["proteins"]
    contigs_df.reset_index(inplace=True)

    #print("******* alla fine di tutto")
    #print(gene2genomeTMP.head())
    #gene2genome=gene2genomeTMP.copy(deep=True)
    return gene2genomeTMP, clusters_df, profiles_df, contigs_df


def load_mcl_clusters(fi):

    """
    Load given clusters file
    
    Args:
        fi (str): path to clusters file
        proteins_df (dataframe): A dataframe giving the protein and its contig.
    Returns: 
        tuple: dataframe proteins and dataframe clusters
    """
    
    # Read MCL
    with open(fi) as f:
        c = [line.rstrip("\n").split("\t") for line in f]

    c = [x for x in c if len(c) > 1]
    nb_clusters = len(c)
    formatter = "PC_{{:>0{}}}".format(int(round(np.log10(nb_clusters))+1))
    name = [formatter.format(str(i)) for i in range(nb_clusters)]
    size = [len(i) for i in c]
    clusters_df = pd.DataFrame({"size": size, "pc_id": name}).set_index("pc_id")
    
    return clusters_df, name, c


def load_one_clusters(fi):
    """
    Load given clusters file

    Args:
        fi (str): path to clusters file
        proteins_df (dataframe): A dataframe giving the protein and its contig.
    Returns:
        tuple: dataframe proteins and dataframe clusters
    """

    fi_clusters_df = pd.read_csv(fi, header=0)

    c = [line.rstrip().split() for line in fi_clusters_df['Members']]
    c = [x for x in c if len(c) > 1]  #
    nb_clusters = len(c)
    formatter = "PC_{{:>0{}}}".format(int(round(np.log10(nb_clusters)) + 1))
    name = [formatter.format(str(i)) for i in range(nb_clusters)]
    size = [len(i) for i in c]

    clusters_df = pd.DataFrame({"size": size, "pc_id": name}).set_index("pc_id")

    return clusters_df, name, c
