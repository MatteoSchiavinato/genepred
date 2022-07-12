#!/usr/bin/env python3

import sys
from Bio import SeqIO as seqio
from math import floor
import argparse as ap

# parser

p = ap.ArgumentParser()
p.add_argument("--input-fasta", required=True)
p.add_argument("--num-out-files", default=16, type=int)
args = p.parse_args()


# functions

def count_sequences(fasta):

    # count sequences
    Sequences = seqio.parse(fasta, "fasta")
    num_seqs = int(len(list(Sequences)))
    sys.stderr.write("Read {0} sequences\n".format(num_seqs))

    return num_seqs



if __name__ == "__main__":

    # read input file
    infile = args.input_fasta
    num_seqs = count_sequences(infile)

    if num_seqs >= args.num_out_files:
        num_out_files = args.num_out_files
        seqs_per_file = floor(num_seqs / num_out_files)
    else:
        num_out_files = num_seqs
        seqs_per_file = 1

    # total processed sequences
    j = 0
    # sequences put in a file
    i = 0
    # file counter
    k = 1

    OUTPUT = open("{0}.{1}.fa".format(infile, k), "w")
    for record in seqio.parse(infile, "fasta"):
        # if the global counter is still not at the end of the iterator
        if ((i < seqs_per_file) and (j <= num_seqs)):
            # simply write to output
            OUTPUT.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")
            i += 1

        # if we've come to the max num sequences
        elif ((i >= seqs_per_file) and (j <= num_seqs)):
            if (k < num_out_files):
                # close previous file
                OUTPUT.close()
                # open new file
                k += 1
                OUTPUT = open("{0}.{1}.fa".format(infile, k), "w")
                # write to new output file
                OUTPUT.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")
                # update counters
                i = 1

            # if this is the last file we put all that remains here
            else:
                OUTPUT.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")
                i += 1

        # global number of processed sequences
        j += 1
        if (j % 1000 == 0):
            # print number to monitor status
            sys.stderr.write("Processed {0} sequences\n".format(j))

    # final else
    else:
        sys.stderr.write("Processed {0} sequences\n".format(j))
        OUTPUT.close()
