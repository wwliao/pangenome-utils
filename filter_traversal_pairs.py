#!/usr/bin/env python3
import argparse
from os.path import basename, splitext

parser = argparse.ArgumentParser()
parser.add_argument("pairs", help="Traversal pairs in BEDPE format")
args = parser.parse_args()

root, ext = splitext(basename(args.pairs))
with open(args.pairs) as infile, open(f"{root}.filtered{ext}", "w") as outfile:
    for i, line in enumerate(infile):
        if i == 0:
            outfile.write(line)
        else:
            filter = line.strip().split("\t")[10]
            if filter == ".":
                outfile.write(line)
