#!/usr/bin/env python3
import glob
import argparse
import statistics
from Bio import SeqIO
from collections import defaultdict


# Author: Peter Mulhair
# Date: 29/03/2022

'''
Script to calculate GC content per
chromosome/scaffold in a genome
with a given cutoff size.
'''

parse = argparse.ArgumentParser()

parse.add_argument("-p", "--path",type=str, help="path to genome fasta files",required=True)
parse.add_argument("-c", "--cutoff",type=int, help="scaffold size cutoff",required=True)
parse.add_argument("-o", "--outfile",type=str, help="name of output file",required=True)
parse.add_argument("-i", "--info",type=str, help="tsv file with species phylum information",required=False)

args = parse.parse_args()

#If species-clade info file provided, store data
if args.info:
    sp_names = {}
    sp_group = {}
    with open(args.info) as f:
        for line in f:
            lines = line.split('\t')
            species = lines[0]
            group = lines[2].strip()
            sp_short= lines[1]
            sp_names[sp_short] = species
            sp_group[sp_short] = group

genome_list = []            
for fasta in glob.glob(args.path + '*'):
    if fasta.endswith(('.fa', '.fasta', '.fas', '.fna')):
            genome_list.append(fasta)
print('Calculating GC content for', len(genome_list), 'genomes...\n')
            
outF = open(args.outfile, 'w')
for fasta in glob.glob(args.path + '*'):
    if fasta.endswith(('.fa', '.fasta', '.fas', '.fna')):
        genome = fasta.split('/')[-1]
        sp_name = genome.split('.')[1].split('_')[1]
        if args.info:
            species = sp_names[sp_name]
            group = sp_group[sp_name]
        
        chr_GC = {}
        print(sp_name)
        with open(fasta) as f:
            for record in SeqIO.parse(f, 'fasta'):
                header = record.description
                seq = str(record.seq)
                seq_len = len(seq)
                GC_count = seq.count('G') + seq.count('g') + seq.count('C') + seq.count('c')
                GC_content = (GC_count/seq_len)*100#Measure GC content by dividing number of Gs and Cs by the scaffold length
                chr_GC[header] = GC_content
                if seq_len >= args.cutoff:
                    if args.info:
                        outF.write(species + '\t' + group + '\t' + str(seq_len) + '\t' + str(GC_content) + '\n')
                    else:
                        outF.write(sp_name + '\t' + str(seq_len) + '\t' + str(GC_content) + '\n')
                    
outF.close()
