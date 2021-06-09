#!/usr/bin/env python3
import re
import argparse
from collections import defaultdict
from os.path import basename, splitext

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

def cvt2stable(pri, traversal):
    if traversal != "*":
        a = []
        for m in re.finditer("([><])([^\s><]+)", traversal):
            if m:
                dir = m[1]
                seg = m[2]
                if seg not in pri:
                    if dir == ">":
                        seq = nonpri_seq[seg]
                    else:
                        seq = revcomp(nonpri_seq[seg])
                    add_new = True
                    if len(a) > 0:
                        b = a[-1]
                        if len(b) == 2:
                            b[1] += seq
                            add_new = False
                    if add_new:
                        a.append(["+", seq])
                else:
                    h = pri[seg]
                    add_new = True
                    if len(a) > 0:
                        b = a[-1]
                        if len(b) == 4:
                            if b[0] == dir and h[0] == b[1]:
                                if b[0] == ">":
                                    if h[1] == b[3]:
                                        b[3] = h[2]
                                        add_new = False
                                else:
                                    if h[2] == b[2]:
                                        b[2] = h[1]
                                        add_new = False

                    if add_new:
                        a.append([dir, h[0], h[1], h[2]])
        stable = []
        for i in range(len(a)):
            if len(a[i]) != 4:
                stable.append("".join(a[i]))
            else:
                stable.append(a[i][0] + a[i][1] + ':' + str(a[i][2]) + '-' + str(a[i][3]))
        return "".join(stable)
    else:
        return "*"

parser = argparse.ArgumentParser()
parser.add_argument("gfa")
parser.add_argument("dipbedpe")
args = parser.parse_args()

seg_pattern = re.compile("^S\t(\S+)\t(\S+)(\t.*)")
tag_pattern = re.compile("\t(SN|SO|SR):[Zi]:(\S+)")
pri = {}
pri_len = defaultdict(int)
nonpri_seq = {}
with open(args.gfa) as infile:
    for line in infile:
        if line.startswith("S"):
            m = seg_pattern.search(line)
            if m:
                seg = m[1]
                if m[2] == "*":
                    length = 0
                else:
                    length = len(m[2])
                tags = m[3]
                for m in tag_pattern.finditer(tags):
                    if m[1] == "SN":
                        sn = m[2]
                    elif m[1] == "SO":
                        so = int(m[2])
                    elif m[1] == "SR":
                        sr = int(m[2])
                if sr == 0:
                    pri[seg] = [sn, so, so + length]

                    if pri_len[sn] < so + length:
                        pri_len[sn] = so + length
            else:
                cols = line.strip().split("\t")
                seg = cols[1]
                seq = cols[2]
                nonpri_seq[seg] = seq

root, ext = splitext(basename(args.dipbedpe))
with open(args.dipbedpe) as infile, open(f"{root}.dip_stable{ext}", "w") as outfile:
    for i, line in enumerate(infile):
        cols = line.strip().split("\t")
        if i == 0:
            cols.insert(9, "Stable_Ref")
            cols.insert(10, "Stable_Alt")
            outfile.write("\t".join(cols) + "\n")
        else:
            outfile.write("\t".join(cols[:9]))
            ref = cols[7]
            stable_ref = cvt2stable(pri, ref)
            outfile.write(f"\t{stable_ref}")
            stable_alts = []
            alts = cols[8].split(",")
            for alt in alts:
                stable_alts.append(cvt2stable(pri, alt))
            stable_alt_str = ",".join(stable_alts)
            outfile.write(f"\t{stable_alt_str}\t")
            outfile.write("\t".join(cols[9:]) + "\n")
