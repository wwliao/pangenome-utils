#!/usr/bin/env python3
import re
import argparse
from os.path import splitext, basename
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--striptags", action="store_true")
parser.add_argument("gfafile")
args = parser.parse_args()

root, ext = splitext(basename(args.gfafile))
counter = Counter()
with open(args.gfafile) as infile, open(f"{root}.w2p{ext}", "w") as outfile:
    for line in infile:
        if line.startswith("W"):
            cols = line.strip().split("\t")
            contig = f"{cols[1]}#{cols[2]}#{cols[3]}"
            counter.update([contig])
            count = counter[contig]
            pathname = f"{contig}#{count}"
            walk = cols[6]
            path = []
            for match in re.finditer("([><])(\d+)", walk):
                node = match.group(2)
                if match.group(1) == ">":
                    direct = "+"
                else:
                    direct = "-"
                path.append(f"{node}{direct}")
            pathstr = ",".join(path)
            outfile.write(f"P\t{pathname}\t{pathstr}\t*\n")
        elif args.striptags and line.startswith("S"):
            segment = "\t".join(line.strip().split("\t")[:3]) + "\n"
            outfile.write(segment)
        elif args.striptags and line.startswith("L"):
            link = "\t".join(line.strip().split("\t")[:6]) + "\n"
            outfile.write(link)
        else:
            outfile.write(line)
