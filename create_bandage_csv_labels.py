#!/usr/bin/env python3
# CSV format: https://github.com/rrwick/Bandage/wiki/CSV-labels
import re
import argparse
from os.path import basename, splitext

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rank", type=int, default=0)
parser.add_argument("-c", "--color", default="#4CB391")
parser.add_argument("gfafile")
args = parser.parse_args()

prefix = splitext(basename(args.gfafile))[0]
with open(args.gfafile) as infile, open(f"{prefix}.labels.csv", "w") as outfile:
    outfile.write("Node,Color\n")
    pattern = re.compile("^S\t(\S+)\t(\S+)(\t.*)*")
    for line in infile:
        color = "#BEBEBE"
        m = pattern.search(line)
        if m:
            node = m[1]
            tags = m[3]
            if tags:
                if re.search(f"\tSR:i:{args.rank}", tags):
                    color = args.color
            outfile.write(f"{node},{color}\n")
