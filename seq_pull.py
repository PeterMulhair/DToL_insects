#!/usr/bin/env python3
import argparse
from Bio import SeqIO

parse = argparse.ArgumentParser()

parse.add_argument("-g", "--genome", type=str, help="genome fasta file to pull sequence from")
parse.add_argument("-c", "--contig", type=str, help="name of contig of interest in genome")
parse.add_argument("-r", "--ranges", type=str, help="range of sequnce to pull out")
parse.add_argument("-o", "--out", type=str, help="name of output file")

args = parse.parse_args()

seq_start = args.ranges.split(',')[0]
seq_end = args.ranges.split(',')[1]
query_range = range(int(seq_start) - 1, int(seq_end))

with open(args.genome) as f, open(args.out, 'a+') as outF:
    for record in SeqIO.parse(f, 'fasta'):
        ID = record.id
        if ID == args.contig:
            seq = str(record.seq)
            query_seq = []
            for i in query_range:
                try:
                    nuc = seq[i]
                    query_seq.append(nuc)
                except:
                    sys.exit('Error: query range not in genome. Double check your input range')
            query_seq = ''.join(query_seq)
                    
            outF.write('>' + ID + ' ' + seq_start + '_' + seq_end + '\n' + query_seq + '\n')
