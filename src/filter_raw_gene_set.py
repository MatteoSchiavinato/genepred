#!/usr/bin/env python3

# modules
import argparse as ap
import sys

# parser
p = ap.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--remove-genes")
p.add_argument("--min-evidence", type=float, default=1.0)
args = p.parse_args()

# functions
def get_id_and_attrs(line):

    lst = line.rstrip("\n\r\b").split("\t")
    features = {x.split("=")[0]:x.split("=")[1] for x in lst[8].split(";")}

    try:
        id = features["ID"]
    except KeyError:
        id = "NA"

    try:
        name = features["Name"]
    except KeyError:
        name = "NA"

    try:
        parent = features["Parent"]
    except KeyError:
        parent = "NA"

    return (id, name, parent)


def get_evidence_for_transcripts(Lines, args):

    Genes_w_evidence = []

    for line in Lines:
        lst = line.rstrip("\b\r\n").split("\t")
        feat = lst[2]

        if feat == "transcript":
            id, name, parent = get_id_and_attrs(line)
            evidence = float({ x.split("=")[0]:x.split("=")[1] for x in lst[8].split(";")}["Note"].replace("evidence_", "").replace("%", ""))

            if evidence >= args.min_evidence:
                Genes_w_evidence.append(parent)
                Genes_w_evidence.append(id)

    return Genes_w_evidence


def find_genes_with_valid_transcripts(Lines):

    Genes_with_transcript = []

    for line in Lines:
        lst = line.rstrip("\b\r\n").split("\t")
        feat = lst[2]

        if feat == "transcript":
            id, name, parent = get_id_and_attrs(line)
            Genes_with_transcript.append(parent)

    return Genes_with_transcript


if __name__ == "__main__":

    # read input file
    IN = open(args.input, "r")
    Lines = [line for line in IN]
    IN.close()

    # read genes to eliminate
    if args.remove_genes:
        IN = open(args.remove_genes, "r")
        Genes_to_remove = [line.rstrip("\b\r\n") for line in IN]
        IN.close()

        # remove genes we don't want
        Filt_lines = []
        for line in Lines:
            id, name, parent = get_id_and_attrs(line)
            if (id not in Genes_to_remove) and (name not in Genes_to_remove) and (parent not in Genes_to_remove):
                Filt_lines.append(line)

        Lines = Filt_lines

    else:
        pass

    # get evidence for each gene
    Genes_w_evidence = get_evidence_for_transcripts(Lines, args)

    # remove transcripts with low evidence
    Filt_lines = []
    for line in Lines:
        id, name, parent = get_id_and_attrs(line)

        if ((id in Genes_w_evidence) \
        or (name in Genes_w_evidence) \
        or (parent in Genes_w_evidence)):
            Filt_lines.append(line)

    Lines = Filt_lines

    # remove genes with no transcripts after evidence removal
    Genes_with_transcript = find_genes_with_valid_transcripts(Lines)

    for line in Lines:
        id, name, parent = get_id_and_attrs(line)
        feat = line.rstrip("\b\r\n").split("\t")[2]

        if feat == "gene":
            gene = id
        elif feat == "transcript":
            gene = parent
        else:
            if id != "NA":
                gene = id.split(".")[0]
            else:
                if name != "NA":
                    gene = name.split(".")[0]
                else:
                    if parent != "NA":
                        gene = parent.split(".")[0]
                    else:
                        sys.stderr.write("This line is weird:\n")
                        sys.stderr.write(line)
                        sys.exit()

        if gene in Genes_with_transcript:
            sys.stdout.write(line)
