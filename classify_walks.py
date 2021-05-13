#!/usr/bin/env python3
import re
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("gfafile")
args = parser.parse_args()

ref_nodes = {}
sample_nodes = defaultdict(lambda: defaultdict(set))
sample_walks = defaultdict(list)
with open(args.gfafile) as infile:
    for line in infile:
        if line.startswith("P"):
            cols = line.strip().split("\t")
            chrom = cols[1]
            path = cols[2]
            ref_nodes[chrom] = set(re.findall("(\d+)[+-],*", path))
        if line.startswith("W"):
            cols = line.strip().split("\t")
            sample = cols[1]
            if sample == args.sample:
                hap = int(cols[2])
                contig = cols[3]
                walk = cols[6]
                sample_nodes[hap][contig].update(re.findall("[><](\d+)", walk))
                sample_walks[contig].append(walk)

with open(f"{args.sample}.0.walks.txt", "w") as nreffile:
    nreffile.write("Chrom\tContig\tWalk\n")
    for hap in sample_nodes:
        with open(f"{args.sample}.{hap}.walks.txt", "w") as outfile:
            outfile.write("Chrom\tContig\tWalk\n")
            for contig in sample_nodes[hap]:
                refcount = 0
                nonref = True
                for chrom in ref_nodes:
                    if len(sample_nodes[hap][contig] & ref_nodes[chrom]) > 0:
                        nonref = False
                        for i, walk in enumerate(sample_walks[contig]):
                            outfile.write(f"{chrom}\t{contig}_{i}\t{walk}\n")
                if nonref:
                    for i, walk in enumerate(sample_walks[contig]):
                        nreffile.write(f"*\t{contig}_{i}\t{walk}\n")
