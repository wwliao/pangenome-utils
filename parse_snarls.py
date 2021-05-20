#!/usr/bin/env python3
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", required=True)
parser.add_argument("ref_nodes")
parser.add_argument("snarls")
args = parser.parse_args()

ref = {}
with open(args.ref_nodes) as infile:
    for line in infile:
        cols = line.strip().split("\t")
        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        node = cols[3]
        ref[node] = (chrom, start, end)

with open(args.snarls) as infile, open(args.output, "w") as outfile, open("filterout.bedpe", "w") as filterout:
    for line in infile:
        snarl = json.loads(line)
        if "parent" not in snarl and "type" in snarl:
            if snarl["type"] == 1:
                node1 = snarl["start"]["node_id"]
                node2 = snarl["end"]["node_id"]
                if node1 in ref and node2 in ref:
                    if "backward" not in snarl["start"] and "backward" not in snarl["end"]:
                        source = node1
                        sink = node2
                        snarl_id = f">{source}>{sink}"
                    elif "backward" in snarl["start"] and "backward" in snarl["end"]:
                        source = node2
                        sink = node1
                        snarl_id = f">{source}>{sink}"
                    else:
                        source = node1
                        sink = node2
                        if "backward" in snarl["start"]:
                            snarl_id = f"<{source}>{sink}"
                        else:
                            snarl_id = f">{source}<{sink}"

                    chrom1 = ref[source][0]
                    start1 = ref[source][1]
                    end1 = ref[source][2]
                    chrom2 = ref[sink][0]
                    start2 = ref[sink][1]
                    end2 = ref[sink][2]
                    if chrom1 == chrom2 and start2 >= end1:
                        outfile.write(f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{snarl_id}\n")
                    else:
                        filterout.write(f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\t{snarl_id}\n")
