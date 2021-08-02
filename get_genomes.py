#!/usr/bin/env python3

import sys
import glob
import argparse
from Bio import SeqIO
from collections import defaultdict
from subprocess import call as unix
from joblib import Parallel, delayed

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage: python3 get_genomes.py --ena <ena tsv file> --group <Name of insect orders>

'''
Script to download available insect DToL
genomes from ncbi. Input can be all insect
genomes or from select insect orders. A 
tsv file from ENA is required as input.
'''

parse = argparse.ArgumentParser()

parse.add_argument("-e", "--ena",type=str, help="ena tsv file with genome info",required=True)
parse.add_argument("-g", "--group", nargs="+", help="list of insect orders to download - eg. Insect, Odonata, Ephemeroptera, Coleoptera, Hymenoptera, Lepidoptera...",required=True)
parse.add_argument("-i", "--info", help="flag to output information on genome data - eg. genome size, chromosome count", action="store_true")

args = parse.parse_args()

#Dictionary of insect orders and DToL abbreviations for them
insect_order_dict = {'Archaeognatha': 'ia', 'Blattodea': 'ib', 'Coleoptera': 'ic', 'Dermaptera': 'ip', 'Diptera': 'id', 'Embioptera': 'ie', 'Ephemeroptera': 'ie', 'Hemiptera': 'ih', 'Hymenoptera': 'iy', 'Lepidoptera': 'il', 'Mantodea': 'im', 'Mecoptera': 'im', 'Megaloptera': 'im', 'Neuroptera': 'in', 'Odonata': 'io', 'Orthoptera': 'io', 'Phasmatodea': 'ip', 'Phthiraptera': 'ip', 'Plecoptera': 'ip', 'Poduromorpha': 'ip', 'Psocoptera': 'ip', 'Raphidioptera': 'ir', 'Siphonaptera': 'is', 'Strepsiptera': 'is', 'Thysanoptera': 'it', 'Trichoptera': 'it', 'Zygentoma': 'iz'}
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

#Parse ena file to create dictionary of species name to GCA assembly ID
genome_dict = {}
genome_info = defaultdict(list)
with open(args.ena) as f:
    next(f)
    for line in f:
        lines = line.split('\t')
        sp = lines[1].strip()
        if 'alternate' not in sp:
            GCA = lines[0].strip('"')
            GCAid = GCA.split('.')[0]
            if GCAid not in GCA_list:
                spID = sp.split(' ')[0].strip('"')
                sp_name = sp.split(' ')[-2:]
                sp_name = ' '.join(sp_name)
                sp_name = sp_name.strip('"')
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

print('\n')
print('Downloading', len(genome_dict), 'genome(s)')
print('This may take some time...')

#Function to download genomes using dictionary of species to assembly IDs as input                   
def genome_download(species, genome):
    ID = genome.split('.')[0].split('_')[1]
    ID1 = ID[:3]
    ID2 = ID[3:6]
    ID3 = ID[6:9]
    unix('sudo wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/' + ID1 + '/' + ID2 + '/' + ID3 + '/' + genome + '_' + species + '/' + genome + '_' + species + '_genomic.fna.gz', shell=True)
    
    unix('sudo gzip -d ' + genome + '_' + species + '_genomic.fna.gz', shell=True)
    unix('sudo mv ' + genome + '_' + species + '_genomic.fna ' + genome + '_' + species + '_genomic.fasta', shell=True)

#Run function in parallel to download multiple genomes at once - currently set to download 10 at a time
Parallel(n_jobs=10)(delayed(genome_download)(k, v) for k, v in genome_dict.items())


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
