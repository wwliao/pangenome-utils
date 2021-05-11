#!/usr/bin/env python3
import re
import argparse
from functools import partial
from collections import Counter
from multiprocessing import Pool

def reverse_path(path):
    d = str.maketrans("><", "<>")
    rev_path = []
    for node in re.findall("[><]\d+", path)[::-1]:
        rev_path.append(node.translate(d))
    return "".join(rev_path)

def search_allele(aln, pattern):
    m = pattern.search(aln)
    if m is not None:
        return m[0]
    else:
        rev_aln = reverse_path(aln)
        m = pattern.search(rev_aln)
        if m is not None:
            return m[0]
    return None

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--threads", type=int, default=1)
parser.add_argument("snarls")
parser.add_argument("gfafile")
args = parser.parse_args()

walks = []
with open(args.gfafile) as infile:
    for line in infile:
        if line.startswith("W"):
            walks.append(line.strip().split("\t")[6])

with open(f"{args.snarls}.alleles.txt", "w") as outfile:
    outfile.write("Source_Node\tSink_Node\tLevel\tParent_Snarl\tAllele\tWalk_Depth\n")
    with open(args.snarls) as infile:
        for i, line in enumerate(infile):
            if i > 0:
                line = line.strip()
                outfile.write(line)
                source, sink = line.split("\t")[:2]
                bubble = (source, sink)
                with Pool(args.threads) as pool:
                    pattern = re.compile(f"{source}(?:[><]\d+)*{sink}")
                    partial_func = partial(search_allele, pattern=pattern)
                    matches = pool.map(partial_func, walks)
                depth = Counter(filter(lambda x: x is not None, matches))
                if depth:
                    alleles, depths = zip(*depth.most_common())
                    allele_str = ",".join(alleles)
                    depth_str = ",".join(map(str, depths))
                    outfile.write(f"\t{allele_str}\t{depth_str}\n")
                else:
                    outfile.write("\t.\t0\n")
