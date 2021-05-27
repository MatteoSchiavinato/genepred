#!/usr/bin/env python3.6

# file columns
# this is a custom output of the "run_blast.sh" script
# this does NOT correspond to the regular blast output

# 1.  qseqid
# 2.  sseqid
# 3.  pident
# 4.  qcovs
# 5.  evalue
# 6.  bitscore
# 7.  staxid
# 8.  ssciname
# 9.  scomname
# 10. sskingdom
# 11. length
# 12. nident
# 13. mismatch
# 14. gapopen
# 15. gaps
# 16. qlen
# 17. qstart
# 18. qend
# 19. qseq

import sys
import pandas as pd
import argparse as ap

# parser
p = ap.ArgumentParser()
p.add_argument("--input-file", help="Output of BLAST format 6 with following columns: qseqid sseqid pident qcovs evalue bitscore staxid ssciname scomname sskingdom length nident mismatch gapopen gaps qlen qstart qend qseq [default: stdin]")
p.add_argument("--max-evalue", default=0.001, type=float, help="Maximum e-value threshold, only hits with at most that evalue will be retained [default: 0.001]")
p.add_argument("--min-pident", default=0, type=float, help="Minimum % of sequence identity from the blast search [Default: off]")
p.add_argument("--min-qcovs", default=0, type=float, help="Minimum query coverage (%) from the blast hit [Default: off]")
p.add_argument("--min-bitscore", default=0, type=float, help="Minimum bistcore from the blast hit [Default: off]")
p.add_argument("--max-mismatches", default=float("inf"), type=float, help="Maximum number of tolerated mismatches from the blast hit [Default: off]")
p.add_argument("--max-gaps", default=float("inf"), type=float, help="Maximum number of tolerated gap openings from the blast hit [Default: off]")
args = p.parse_args()


# input and output
if args.input_file:
    INPUT = open(sys.argv[1], "r")
else:
    INPUT = sys.stdin

OUTPUT = sys.argv[1] + "." + "filt"

# save all results as a dictionary
Lines = {}
for line in INPUT:
    lst = line.rstrip("\n\r\b").split("\t")
    gene = lst[0]
    try:
        Lines[gene].append(lst)
    except KeyError:
        Lines[gene] = [lst]


# treat each dictionary as a pandas dataframe
# sort by the columns needed (evalue 5(#4) > pident 3(#2) > qcovs 4(#3) > bitscore 6(#5))
# take the first row after sorting
Best_hits = pd.DataFrame(columns = ["Query", "Subject", "Pident", "Qcovs", "Evalue",
                                    "Bitscore", "Taxid", "Sciname", "Comname", "Kingdom",
                                    "Length", "Nident", "Mismatch", "Gapopen", "Gaps",
                                    "Qlen", "Qstart", "Qend", "Qseq", "Stitle" ])

for gene in Lines.keys():
    df = pd.DataFrame(Lines[gene])
    df = df.drop_duplicates()
    # filter hits based on criteria
    df = df[df.Evalue <= args.max_evalue , ]
    df = df[df.Pident >= args.min_pident , ]
    df = df[df.Qcovs >= args.min_qcovs , ]
    df = df[df.Bitscore >= args.min_bitscore , ]
    df = df[df.Gaps <= args.max_gaps , ]
    df = df[df.Mismatch <= args.max_mismatches , ]
    # sort and select best
    df = df.sort_index(axis=0, level=[4,2,3,5])
    df.columns = ["Query", "Subject", "Pident", "Qcovs", "Evalue",
                  "Bitscore", "Taxid", "Sciname", "Comname", "Kingdom",
                  "Length", "Nident", "Mismatch", "Gapopen", "Gaps",
                  "Qlen", "Qstart", "Qend", "Qseq", "Stitle" ]
    best_hit = pd.Series(df.iloc[0,])
    Best_hits = Best_hits.append(df.iloc[[0], : ])


Best_hits.to_csv(OUTPUT, sep="\t", index=False)
