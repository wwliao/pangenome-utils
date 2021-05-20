#!/usr/bin/env python3
import re
import argparse
import json
from collections import defaultdict, Counter

def reverse_path(path):
    d = str.maketrans("><", "<>")
    rev_path = []
    for node in re.findall("[><]\d+", path)[::-1]:
        rev_path.append(node.translate(d))
    return "".join(rev_path)

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", required=True)
parser.add_argument("snarls", help="Snarls in BEDPE format")
parser.add_argument("traversals", help="Traversals in JSON format")
args = parser.parse_args()

snarls = set()
with open(args.snarls) as infile:
    for line in infile:
        cols = line.strip().split("\t")
        snarl = cols[6]
        snarls.add(snarl)

traversals = defaultdict(list)
with open(args.traversals) as infile:
    for line in infile:
        traversal = []
        for step in json.loads(line)["visit"]:
            if "backward" not in step:
                traversal.append(f">{step['node_id']}")
            else:
                traversal.append(f"<{step['node_id']}")
        source = traversal[0]
        sink = traversal[-1]
        snarl = source + sink
        rev_snarl = reverse_path(snarl)
        if snarl in snarls:
            traversals[snarl].append("".join(traversal))
        elif rev_snarl in snarls:
            rev_traversal = reverse_path("".join(traversal))
            traversals[rev_snarl].append(rev_traversal)

with open(args.snarls) as infile, open(args.output, "w") as outfile:
    outfile.write("Chrom1\tStart1\tEnd1\tChrom2\tStart2\tEnd2\tSnarl\tTraversal\tContig_Coverage\n")
    for line in infile:
        line = line.strip()
        snarl = line.split("\t")[6]
        outfile.write(line)
        if len(traversals[snarl]) > 0:
            counter = Counter(traversals[snarl])
            alleles = []
            counts = []
            for allele, count in counter.most_common():
                alleles.append(allele)
                counts.append(str(count))
            allele_str = ",".join(alleles)
            count_str = ",".join(counts)
            outfile.write(f"\t{allele_str}\t{count_str}\n")
        else:
            outfile.write(f"\t.\t0\n")
