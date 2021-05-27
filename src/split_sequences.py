#!/usr/bin/env python3.6

import sys
from Bio import SeqIO as seqio

# count sequences
infile = sys.argv[1]
Sequences = seqio.parse(infile, "fasta")
num_seqs = int(len(list(Sequences)))
sys.stderr.write("Read {0} sequences\n".format(num_seqs))

# sequences in a file
i = 1
# number of files
k = 1
# total processed sequences
j = 0

OUTPUT = open("{0}.{1}.fa".format(infile, k), "w")
for record in seqio.parse(infile, "fasta"):
    # if the global counter is still not at the end of the iterator
    if (j <= num_seqs):
        # if we can still fill up the containing split file
        if (i < 1000):
            # simply write to output
            OUTPUT.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")
            i += 1
        # if we've come to 1000 sequences
        elif (i == 1000):
            # write to output and close file
            OUTPUT.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")
            OUTPUT.close()
            # update counters
            i = 1
            k += 1
            # reopen new output file
            OUTPUT = open("{0}.{1}.fa".format(infile, k), "w")
        else:
            sys.stderr.write("ERROR: we don't want this to happen (i > 1000)\n")
    else:
        sys.stderr.write("ERROR: we don't want this to happen (j > len(Sequences))\n")

    # global number of processed sequences
    j += 1
    if (j % 1000 == 0):
        # print number to monitor status
        sys.stderr.write("Processed {0} sequences\n".format(j))

# final else
else:
    sys.stderr.write("Processed {0} sequences\n".format(j))
    OUTPUT.close()
