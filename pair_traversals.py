#!/usr/bin/env python3
import argparse

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
                outfile.write("\tRef\tPat\tMat\tFilter\tContig_Coverage\tAverage_Edge_Depth\n")
            else:
                filter = []
                ref_traversal = ref_line.split("\t")[7]
                pat_cols = pat_line.split("\t")
                pat_traversal = pat_cols[7]
                pat_cov = pat_cols[8]
                pat_depth = pat_cols[9]
                mat_cols = mat_line.split("\t")
                mat_traversal = mat_cols[7]
                mat_cov = mat_cols[8]
                mat_depth = mat_cols[9]
                outfile.write("\t".join(ref_line.split("\t")[:8]))
                outfile.write(f"\t{pat_traversal}\t{mat_traversal}")

                if ref_traversal == ".":
                    filter.append("GAP0")

                if pat_traversal == ".":
                    filter.append("GAP1")
                elif len(pat_traversal.split(",")) > 1:
                    filter.append("HET1")

                if mat_traversal == ".":
                    filter.append("GAP2")
                elif len(mat_traversal.split(",")) > 1:
                    filter.append("HET2")
                if len(filter) > 0:
                    filter_str = ";".join(filter)
                else:
                    filter_str = "."
                outfile.write(f"\t{filter_str}\t{pat_cov}:{mat_cov}\t{pat_depth}:{mat_depth}\n")
