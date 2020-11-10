#!/usr/bin/env python

import sys
import csv
from pyparsing import Literal, SkipTo
import argparse

parser = argparse.ArgumentParser(
    description='This script is designed to prepare a gene2contig file, appropriate for use with vContact-PCs '
                '(i.e. proteins.csv) using a variety of different input formats. It currently accepts output from '
                'VIRSorter, Prodigal, MetaGeneMark and NCBI\'s coding features fasta format. While this generates a'
                'gene2contig file, there is not sufficient information in any of the above format\'s headers to append'
                'annotation information.',
    formatter_class=argparse.RawTextHelpFormatter)

options = parser.add_argument_group('Options')
options.add_argument('-p', '--proteins', dest='proteins', metavar='FILENAME', default='VIRSorter_prots.fasta')
options.add_argument('-o', '--output', dest='output', metavar='FILENAME', default='proteins.csv')
options.add_argument('-s', '--source-type', dest='source_type',
                     choices=['VIRSorter', 'Prodigal-coords', 'Prodigal-FAA', 'MetaGeneMark', 'NCBICodingSequence',
                              'NCBIFasta'],
                     required=True,
                     help='Select one of the options as an input source. MetaGeneMark can be either the nucleotide or '
                          'protein FASTA-formatted output.')
options.add_argument('-k', '--keep-descriptions', dest='keep_description', action='store_true',
                     help='This will enable taking the full description of sequences during MetaGeneMark parsing.')
options.add_argument('-c', '--compatibility', dest='compatibility', action='store_true',
                     help='Adds compatibility for vContact1 headers.')

results = parser.parse_args()


def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)


try:
    from Bio import SeqIO
except ImportError:
    error("The biopython library was not found.")

gene2contig = {}


def parseVIRSorter(proteins_fh):
    """"""
    relDict = {}
    for protein_records in SeqIO.parse(proteins_fh, 'fasta'):
        proteinID = protein_records.id
        contig = proteinID.split('-gene_')[0]

        if contig not in relDict:
            relDict[contig] = []

        relDict[contig].append(proteinID)

    return relDict


def parseProdigalCoords(proteins_fh):
    """"""
    relDict = {}
    content = proteins_fh.read()

    start_header = Literal('DEFINITION')
    end_header = Literal('//\n')
    text = start_header + SkipTo(end_header)

    for tokens, start_pos, end_pos in text.scanString(content):  # Tokens[0] is start

        headerDict = dict(ele.split("=") for ele in tokens[1].split('\n')[0].strip().split(';'))

        contig = str(headerDict['seqhdr'].strip('"'))

        if contig not in relDict:
            relDict[contig] = []

        data = tokens[1].strip().split('\n')[1:]
        cdsData = [ele.strip().split('/note="')[-1] for ele in data[2:] if 'CDS' not in ele]

        cdsDict = {}

        for cds in cdsData:
            for datum in cds.split(';'):
                if '"' not in datum:
                    key, value = datum.split('=')
                    cdsDict[key] = value

            if cdsDict['ID'] not in relDict[contig]:
                relDict[contig].append(str(cdsDict['ID']))

    return relDict


def parseMetaGeneMark(proteins_fh):
    """"""
    relDict = {}
    for record in SeqIO.parse(proteins_fh, 'fasta'):

        proteinID = record.id

        if results.keep_description:  # Separate gene from contig with a single split
            contig = str(record.description).split(None, 1)[1]  # 0 is gene
        else:  # Separate ... and don't have anything else back together
            contig = str(record.description).split()[1]

        if contig not in relDict:
            relDict[contig] = []

        relDict[contig].append(proteinID)

    return relDict


def parse_NCBI_coding(proteins_fh):
    """"""
    relDict = {}

    for record in SeqIO.parse(proteins_fh, 'fasta'):
        proteinID = record.id  # lcl|NC_001447.1_prot_NP_040809.1_2
        proteinInfo = record.description.split('protein=')[1].split(']')[0]

        contig_components = proteinID.split('|')[1].split('_')
        contig = '_'.join([contig_components[0], contig_components[1]])

        if contig not in relDict:
            relDict[contig] = []

        relDict[contig].append((proteinID, proteinInfo))

    return relDict


def parseProdigalFAA(proteins_fh):

    relDict = {}

    for record in SeqIO.parse(proteins_fh, 'fasta'):
        # lcl|NC_001447.1_prot_NP_040809.1_2
        proteinID = record.id
        #proteinInfo = record.description.split()[1]

        contig = proteinID.rsplit('_', 1)[0]

        if contig not in relDict:
            relDict[contig] = []

        relDict[contig].append((proteinID))

    return relDict


def parse_NCBI_FAA(proteins_fh):

    relDict = {}

    for record in SeqIO.parse(proteins_fh, 'fasta'):
        proteinID = record.id
        contig = record.description.rsplit('[')[-1].strip(']')
        proteinInfo = record.description.split('[', 1)[0].split(' ', 1)[-1]

        if contig not in relDict:
            relDict[contig] = []

        relDict[contig].append((proteinID, proteinInfo))

    return relDict


if __name__ == "__main__":

    with open(results.proteins, 'r') as proteins_fh:

        if 'VIRSorter' in results.source_type:
            gene2contig = parseVIRSorter(proteins_fh)

        if 'Prodigal-coords' in results.source_type:
            gene2contig = parseProdigalCoords(proteins_fh)

        if 'Prodigal-FAA' in results.source_type:
            gene2contig = parseProdigalFAA(proteins_fh)

        if 'MetaGeneMark' in results.source_type:
            gene2contig = parseMetaGeneMark(proteins_fh)

        if 'NCBICodingSequence' in results.source_type:
            gene2contig = parse_NCBI_coding(proteins_fh)

        if 'NCBIFasta' in results.source_type:
            gene2contig = parse_NCBI_FAA(proteins_fh)

    with open(results.output, 'w') as proteins_fh:

        csvWriter = csv.writer(proteins_fh, delimiter=',', quotechar='"')

        if results.compatibility:
            csvWriter.writerow(['id', 'contig', 'keywords'])
        else:
            csvWriter.writerow(['protein_id', 'contig_id', 'keywords'])

        sources = ['NCBIFasta', 'NCBICodingSequence']
        if any(source in results.source_type for source in sources):
            for contig, genes in gene2contig.items():
                for (gene, info) in genes:
                    csvWriter.writerow([gene, contig, info])
        else:
            for contig, genes in gene2contig.items():
                for gene in genes:
                    csvWriter.writerow([gene, contig, 'None_provided'])
