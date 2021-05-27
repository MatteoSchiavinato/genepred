#!/usr/bin/env python3.6

import argparse as ap
import sys

p = ap.ArgumentParser()
p.add_argument("-i", "--input")
p.add_argument("-o", "--output")
args = p.parse_args()


# reading evidence and transcript names
INPUT = open(args.input, "r")
Evidence = {}
Lines = []

for line in INPUT:
    if line[0:1] != "#":
        lst = line.rstrip("\n\r\b").split("\t")
        Lines.append(lst)
        if lst[2] in ["mRNA", "transcript"]:
            f9 = lst[8].split(";")
            f9dic = { item.split("=")[0]:item.split("=")[1] for item in f9 }
            transcript = f9dic["ID"]
    elif ("# % of transcript supported by hints (any source)" in line):
        evidence = line.rstrip("\n\r\b").split(": ")[1]
        Evidence[transcript] = float(evidence)
        transcript = ""
        evidence = ""

INPUT.close()


# add evidence values to field #9 of the GFF3 file
if args.output:
    OUTPUT = open(args.output, "w")
else:
    OUTPUT = sys.stdout

for lst in Lines:
    if lst[2] in ["mRNA", "transcript"]:
        f9 = lst[8].split(";")
        f9dic = { item.split("=")[0]:item.split("=")[1] for item in f9 }
        transcript = f9dic["ID"]
        evidence = Evidence
        lst[8] = lst[8] + ";" + "Note=evidence_{0}%".format(Evidence[transcript])
        OUTPUT.write("\t".join(lst) + "\n")

    else:
        OUTPUT.write("\t".join(lst) + "\n")
