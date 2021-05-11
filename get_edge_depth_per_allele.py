#!/usr/bin/env python3
import re
import argparse
from statistics import mean
from os.path import basename, splitext

parser = argparse.ArgumentParser()
parser.add_argument("alleles")
parser.add_argument("depths")
args = parser.parse_args()

def reverse_snarl(snarl):
    d = str.maketrans("><", "<>")
    source, sink = snarl
    rev_snarl = (sink.translate(d), source.translate(d))
    return rev_snarl

depths = {}
with open(args.depths) as infile:
    for i, line in enumerate(infile):
        if i > 0:
            cols = line.strip().split("\t")
            source = cols[0]
            sink = cols[1]
            depth = int(cols[2])
            depths[(source, sink)] = depth

root, ext = splitext(basename(args.alleles))
with open(args.alleles) as infile, open(f"{root}.edge_depth{ext}", "w") as outfile:
    pattern = re.compile("([><]\d+)")
    for i, line in enumerate(infile):
        line = line.strip()
        if i == 0:
            header = line
            outfile.write(header)
            outfile.write("\tAverage_Edge_Depth\n")
        else:
            cols = line.split("\t")
            source = cols[0]
            sink = cols[1]
            alleles = cols[4].split(",")
            avg_depths = []
            for allele in alleles:
                nodes = pattern.findall(allele)
                edge_depths = []
                for j in range(len(nodes) - 1):
                    snarl = (nodes[j], nodes[j+1])
                    if snarl in depths:
                        edge_depths.append(depths[snarl])
                    else:
                        edge_depths.append(depths[reverse_snarl(snarl)])
                avg_depths.append(mean(edge_depths))
            avg_depth_str = ",".join(map(lambda d: f"{d:.6f}", avg_depths))
            outfile.write(line)
            outfile.write(f"\t{avg_depth_str}\n")
