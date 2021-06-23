#!/usr/bin/env python3
import re
import argparse
from os.path import splitext, basename
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path-prefix", metavar="NAME", required=True, help="All paths beginning with NAME kept as P-lines")
parser.add_argument("gfafile")
args = parser.parse_args()

root, ext = splitext(basename(args.gfafile))
node_size = {}
path_size = defaultdict(int)
with open(args.gfafile) as infile, open(f"{root}.p2w{ext}", "w") as outfile:
    for line in infile:
        if line.startswith("S"):
            outfile.write(line)
            cols = line.strip().split("\t")
            node = cols[1]
            size = len(cols[2])
            node_size[node] = size
        if line.startswith("P"):
            cols = line.strip().split("\t")
            pid = cols[1]
            if pid.startswith(args.path_prefix):
                outfile.write(line)
            else:
                if len(pid.split("#")) == 2:
                    sample, contig = pid.split("#")
                    sample = sample.upper()
                    hap = "0"
                elif len(pid.split("#")) == 3:
                    sample, hap, contig = pid.split("#")
                else:
                    print(f"{pid}: Path name format is not correct!")
                    break

                start = path_size[pid]
                wline = ""
                for m in re.finditer("(\d+)([+-])", cols[2]):
                    if m:
                        node = m[1]
                        strand = m[2]
                        if strand == "+":
                            dir = ">"
                        else:
                            dir = "<"
                        wline += f"{dir}{node}"
                        path_size[pid] += node_size[node]
                end = path_size[pid]
                outfile.write(f"W\t{sample}\t{hap}\t{contig}\t{start}\t{end}\t{wline}\n")
        else:
            outfile.write(line)
