#!/usr/bin/env python3
import argparse
from Bio import SeqIO

parse = argparse.ArgumentParser()

parse.add_argument("-g", "--genome", type=str, help="genome fasta file to pull sequence from")
parse.add_argument("-c", "--contig", type=str, help="name of contig of interest in genome")
parse.add_argument("-r", "--ranges", action="append", nargs="+", help="range of sequences to pull out eg. 5,10 30,55")
parse.add_argument("-o", "--out", type=str, help="name of output file")

args = parse.parse_args()

with open(args.genome) as f, open(args.out, 'a+') as outF:
    for record in SeqIO.parse(f, 'fasta'):
        ID = record.id
        if ID == args.contig:
            seq = str(record.seq)

            query_seq = []
            for seq_ranges in args.ranges[0]:
                seq_start = seq_ranges.split(',')[0]
                seq_end = seq_ranges.split(',')[1]
                query_range = range(int(seq_start) - 1, int(seq_end))
                for i in query_range:
                    try:
                        nuc = seq[i]
                        query_seq.append(nuc)
                    except:
                        sys.exit('Error: query range not in genome. Double check your input range')
            query_seq = ''.join(query_seq)
            query_start = args.ranges[0][0].split(',')[0]
            query_end = args.ranges[0][1].split(',')[1]
            outF.write('>' + ID + ' ' + query_start + '_' + query_end + '\n' + query_seq + '\n')
