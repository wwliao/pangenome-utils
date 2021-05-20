#!/usr/bin/env python3
import re
import argparse
from os.path import basename, splitext

parser = argparse.ArgumentParser()
parser.add_argument("gfafile")
args = parser.parse_args()

prefix = splitext(basename(args.gfafile))[0]
with open(args.gfafile) as infile, open(f"{prefix}.ref_nodes.bed", "w") as outfile:
    pattern = re.compile("^S\t(\S+)\t(\S+)(\t.*)*")
    for line in infile:
        m = pattern.search(line)
        if m:
            node = m[1]
            seq = m[2]
            if m[3]:
                tags = {}
                for s in re.finditer(f"\t(SR|SN|SO):[iZ]:(\S+)", m[3]):
                    if s[1] in ["SR", "SO"]:
                        tags[s[1]] = int(s[2])
                    else:
                        tags[s[1]] = s[2]
                if tags["SR"] == 0:
                    chrom = tags["SN"]
                    start = tags["SO"]
                    end = start + len(seq)
                    outfile.write(f"{chrom}\t{start}\t{end}\t{node}\n")
