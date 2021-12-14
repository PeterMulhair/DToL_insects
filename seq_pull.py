#!/usr/bin/env python3
import argparse
from Bio import SeqIO

parse = argparse.ArgumentParser()

parse.add_argument("-g", "--genome", type=str, help="genome fasta file to pull sequence from")
parse.add_argument("-c", "--contig", type=str, help="name of contig of interest in genome")
parse.add_argument("-r", "--ranges", action="append", nargs="+", help="range of sequences to pull out eg. 5,10 30,55 (Optional)")
parse.add_argument("-o", "--out", type=str, help="name of output file")

args = parse.parse_args()

with open(args.genome) as f, open(args.out, 'a+') as outF:
    for record in SeqIO.parse(f, 'fasta'):
        ID = record.id
        if ID == args.contig:
            seq = str(record.seq)

            if args.ranges:
                for seq_ranges in args.ranges[0]:
                    query_seq = []
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
                    outF.write('>' + ID + ' ' + seq_start + '_' + seq_end + '\n' + query_seq + '\n')
                
            else:
                outF.write('>' + ID + '\n')
                for seqs in [seq[i:i+80] for i in range(0, len(seq), 80)]:
                    outF.write(str(seqs) + '\n')
