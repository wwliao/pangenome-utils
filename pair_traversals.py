#!/usr/bin/env python3
import re
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", required=True)
parser.add_argument("ref")
parser.add_argument("pat")
parser.add_argument("mat")
args = parser.parse_args()

with open(args.ref) as reffile, open(args.pat) as patfile, open(args.mat) as matfile:
    with open(args.output, "w") as outfile:
        for i, (ref_line, pat_line, mat_line) in enumerate(zip(reffile, patfile, matfile)):
            ref_line = ref_line.strip()
            pat_line = pat_line.strip()
            mat_line = mat_line.strip()
            if i == 0:
                header = "\t".join(ref_line.split("\t")[:7])
                outfile.write(header)
                outfile.write("\tRef\tAlt\tFilter\tGenotype\tEdge_Depth\n")
            else:
                filter = []
                ref_traversal = ref_line.split("\t")[7]
                pat_cols = pat_line.split("\t")
                pat_traversal = pat_cols[7]
                pat_depth = pat_cols[9]
                mat_cols = mat_line.split("\t")
                mat_traversal = mat_cols[7]
                mat_depth = mat_cols[9]

                alleles = [ref_traversal]
                depths = defaultdict(int)
                for traversal, depth in zip(pat_traversal.split(","), pat_depth.split(",")):
                    depths[traversal] = int(depth)
                    if traversal != "." and traversal not in alleles:
                        alleles.append(traversal)
                for traversal, depth in zip(mat_traversal.split(","), mat_depth.split(",")):
                    depths[traversal] = int(depth)
                    if traversal != "." and traversal not in alleles:
                        alleles.append(traversal)

                if pat_traversal == ".":
                    filter.append("GAP1")
                    pat_gt = "."
                elif len(pat_traversal.split(",")) > 1:
                    filter.append("HET1")
                    pat_gt = "."
                else:
                    pat_gt = alleles.index(pat_traversal)

                if mat_traversal == ".":
                    filter.append("GAP2")
                    mat_gt = "."
                elif len(mat_traversal.split(",")) > 1:
                    filter.append("HET2")
                    mat_gt = "."
                else:
                    mat_gt = alleles.index(mat_traversal)

                if len(filter) > 0:
                    filter_str = ";".join(filter)
                else:
                    filter_str = "."

                outfile.write("\t".join(ref_line.split("\t")[:7]))
                ref = re.search("[><]\d+((?:[><]\d+)*)[><]\d+", ref_traversal)[1]
                if ref:
                    outfile.write(f"\t{ref}")
                else:
                    outfile.write("\t*")
                if len(alleles) > 1:
                    alts = []
                    for allele in alleles[1:]:
                        alt = re.search("[><]\d+((?:[><]\d+)*)[><]\d+", allele)[1]
                        if alt:
                            alts.append(alt)
                        else:
                            alts.append("*")
                    alt_str = ",".join(alts)
                    outfile.write(f"\t{alt_str}")
                else:
                    outfile.write("\t.")
                genotype = f"{pat_gt}|{mat_gt}"
                depth_str = ",".join(str(depths[allele]) for allele in alleles)
                outfile.write(f"\t{filter_str}\t{genotype}\t{depth_str}\n")
