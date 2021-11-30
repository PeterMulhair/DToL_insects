#!/usr/bin/env python3

import re
import os.path
import sys
import glob
import argparse
from Bio import SeqIO
import xml.etree.ElementTree as ET
from collections import defaultdict
from subprocess import call as unix
from joblib import Parallel, delayed

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage: python3 get_genomes.py --input <ncbi xml file> --group <Name of insect orders>

'''
Script to download available insect DToL
genomes from ncbi. Input can be all insect
genomes or from select insect orders. A 
tsv file from ENA is required as input.
'''

parse = argparse.ArgumentParser()

parse.add_argument("-i", "--input",type=str, help="ncbi xml file with genome info",required=True)
parse.add_argument("-g", "--group", nargs="+", help="list of insect orders to download - eg. Insect, Odonata, Ephemeroptera, Coleoptera, Hymenoptera, Lepidoptera...",required=True)
parse.add_argument("-d", "--info", help="flag to output information on genome data - eg. genome size, chromosome count", action="store_true")
parse.add_argument("-t", "--threads", help="number of threads to use i.e. number of genomes to download at once, default is 1", required=False)

args = parse.parse_args()

#Dictionary of insect orders and DToL abbreviations for them
insect_order_dict = {'Archaeognatha': 'ia', 'Blattodea': 'ib', 'Coleoptera': 'ic', 'Dermaptera': 'ig', 'Diptera': 'id', 'Embioptera': 'ie', 'Ephemeroptera': 'ie', 'Hemiptera': 'ih', 'Hymenoptera': 'iy', 'Lepidoptera': 'il', 'Mantodea': 'im', 'Mecoptera': 'ij', 'Megaloptera': 'ik', 'Neuroptera': 'in', 'Odonata': 'io', 'Orthoptera': 'iq', 'Phasmatodea': 'ip', 'Phthiraptera': 'ip', 'Plecoptera': 'ip', 'Poduromorpha': 'ip', 'Psocoptera': 'ip', 'Raphidioptera': 'ir', 'Siphonaptera': 'is', 'Strepsiptera': 'is', 'Thysanoptera': 'it', 'Trichoptera': 'ii', 'Zygentoma': 'iz'}
##Collembola instead of Poduromorpha? 

#Check input options for errors
for order_name in args.group:
    if order_name != 'Insecta':
        if order_name not in insect_order_dict:
            print('ERROR: Please check input Order name...')
            print('\n')
            print('Order name should be one of the following:')
            print(list(insect_order_dict.keys()))
            print('Insecta for full download')
            sys.exit()

#Create list of genomes already downloaded to ignore
GCA_list = []
for genomes in glob.glob('*fasta'):
    GCA = genomes.split('.')[0]
    GCA_list.append(GCA)

#Parse ncbi xml file to create dictionary of species name to GCA assembly ID
genome_dict = {}
genome_info = defaultdict(list)
with open(args.input) as f:
    xml = f.read()

ncbi_file = ET.fromstring("<root>" + '\n' + xml + "</root>")
for elem in ncbi_file:

    for subelem in elem.findall('AssemblyName'):
        spID = subelem.text
        
    for subelem in elem.findall('AssemblyAccession'):
        GCA = subelem.text
        GCAid = GCA.split('.')[0]
        if 'GCF' in GCA:
            GCA = 'GCA_' + GCA.split('_')[1]
            GCAid = 'GCA_' + GCAid.split('_')[1]
        
    for subelem in elem.findall('Organism'):
        sp_name = subelem.text
        sp_name = sp_name.split(' (')[0]
        
    if 'alternate' not in spID:
        if GCAid not in GCA_list:
            if args.group[0] == 'Insecta':
                if spID[:1] == 'i':
                    for k, v in insect_order_dict.items():
                        if spID[:2] == v:
                            order_name = k
                    print(sp_name, order_name)
                    genome_dict[spID] = GCA
                    genome_info[order_name].append(sp_name)
            else:
                for order_name in args.group:
                    orderID = insect_order_dict[order_name]
                    if spID[:2] == orderID:
                        print(sp_name, order_name)
                        genome_dict[spID] = GCA
                        genome_info[order_name].append(sp_name)


if len(genome_dict) == 0:
    print('No new genomes to download!')
    sys.exit()
else:
    print('\n')
    print('Downloading', len(genome_dict), 'genome(s)')
    print('This may take some time...')
    print('\n')
    
#Function to download genomes using dictionary of species to assembly IDs as input                   
def genome_download(species, genome):
    ID = genome.split('.')[0].split('_')[1]
    ID1 = ID[:3]
    ID2 = ID[3:6]
    ID3 = ID[6:9]
    unix('sudo wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/' + ID1 + '/' + ID2 + '/' + ID3 + '/' + genome + '_' + species + '/' + genome + '_' + species + '_genomic.fna.gz', shell=True)
    genome_file = genome + '_' + species + '_genomic.fna.gz'
    if os.path.isfile(genome_file):
        unix('sudo gzip -d ' + genome + '_' + species + '_genomic.fna.gz', shell=True)
        unix('sudo mv ' + genome + '_' + species + '_genomic.fna ' + genome + '_' + species + '_genomic.fasta', shell=True)
    else:
        #sys.exit('\nOops, ' + species + 'genome not yet available from ncbi')
        print('Oops, ' + species + 'genome not yet available from ncbi')
    
#Run function in parallel to download multiple genomes at once - use --threads to set how many, default is 1
if args.threads:
    #Get number of threads
    threads = int(args.threads)
    Parallel(n_jobs=threads)(delayed(genome_download)(k, v) for k, v in genome_dict.items())
else:
    Parallel(n_jobs=1)(delayed(genome_download)(k, v) for k, v in genome_dict.items())

#Output useful information on downloaded species if -i/--info flag is given
if args.info:
    print('\n')
    print('Summarizing genome data')
    genome_info_dict = defaultdict(list)
    for fasta in glob.glob('*fasta'):
        sp_short = fasta.split('.')[1].split('_')[1][2:-1]
        order = fasta.split('.')[1].split('_')[1][:2]
        with open(fasta) as f:
            genome_len = 0
            chr_count = 0
            for record in SeqIO.parse(f, 'fasta'):
                header = record.description
                spName = header.split(' ')[1] + ' ' + header.split(' ')[2]
                seq = str(record.seq)
                seq_len = len(seq)
                genome_len+=seq_len
                if 'chromosome' in header:
                    chr_count+=1

        info_list = spName, sp_short, genome_len, chr_count
        genome_info_dict[order].append(info_list)


    with open('genome_summary.tsv', 'a+') as outF:
        for k, v in insect_order_dict.items():
            for orderID, info in genome_info_dict.items():
                if v == orderID:
                    for sp_info in info:
                        outF.write(k + '\t' + sp_info[0] + '\t' + sp_info[1] + '\t' + str(sp_info[2]) + '\t' + str(sp_info[3]) + '\n')


    print('\n')
    print('Download complete. See genome_summary.tsv for useful information about genomes.')
    
else:
    print('\n')
    print('Download complete.')
                        
