#!/usr/bin/env python3
import re
import argparse
from os.path import splitext, basename
from pysam import FastaFile

parser = argparse.ArgumentParser()
parser.add_argument("fasta", help="reference sequences in FASTA format")
parser.add_argument("dipbedpe", help="Diploid stable traversals in BEDPE format")
args = parser.parse_args()

pattern = re.compile("([><][^\s><]+:\d+-\d+)|\+([ACGTacgtRYMKrymkVBHDvbhdWSNwsn]+)")
with FastaFile(args.fasta) as reffile:
    root, ext = splitext(basename(args.dipbedpe))
    with open(args.dipbedpe) as infile, open(f"{root}.classified{ext}", "w") as outfile:
        for i, line in enumerate(infile):
            line = line.strip()
            if i == 0:
                outfile.write(line)
                outfile.write("\tVariant_Type\tVariant_Size\n")
            else:
                cols = line.split("\t")
                if cols[9] != "*":
                    m = re.search(">([^\s><]+):(\d+)-(\d+)", cols[9])
                    chrom = m[1]
                    start = int(m[2])
                    end = int(m[3])
                    refseq = reffile.fetch(chrom, start, end)
                else:
                    refseq = ""
                alts = cols[10].split(",")
                vtypes = []
                vsizes = []
                for j, alt in enumerate(alts):
                    vtypes.append([])
                    vsizes.append([])
                    if alt != "*":
                        steps = []
                        ref_idx = []
                        alt_idx = []
                        for idx, m in enumerate(pattern.finditer(alt)):
                            if m[1]:
                                ref_idx.append(idx)
                                s = re.search("([><])([^\s><]+):(\d+)-(\d+)", m[1])
                                steps.append((s[1], s[2], int(s[3]), int(s[4])))
                            else:
                                alt_idx.append(idx)
                                steps.append(m[2])

                        #if len(alt_idx) <= 1:
                        if len(alt_idx) >= 0:

                            alt_length = 0
                            for idx in alt_idx:
                                alt_length += len(steps[idx])

                            ref_length = 0
                            inv_length = 0
                            isinv = False
                            for idx in ref_idx:
                                if steps[idx][0] != "<":
                                    ref_length += len(reffile.fetch(*steps[idx][1:]))
                                else:
                                    isinv = True
                                    inv_length += len(reffile.fetch(*steps[idx][1:]))

                            if isinv:
                                vtypes[j].append("INV")
                                vsizes[j].append(inv_length)
                                if len(refseq) == ref_length + inv_length + alt_length:
                                    if alt_length > 0:
                                        vtypes[j].append("SNP")
                                        vsizes[j].append(alt_length)
                                else:
                                    if alt_length > 0:
                                        vtypes[j].append("INDEL")
                                        len_diff = (ref_length + inv_length + alt_length) - len(refseq)
                                        vsizes[j].append(len_diff)
                            else:
                                len_diff = (ref_length + alt_length) - len(refseq)
                                if len_diff == 0:
                                    seq = ""
                                    for step in steps:
                                        if type(step) is tuple:
                                            seq += reffile.fetch(*step[1:])
                                        else:
                                            seq += step
                                    if refseq == seq:
                                        vtypes[j].append("REF")
                                        vsizes[j].append(0)
                                    else:
                                        vtypes[j].append("SNP")
                                        vsizes[j].append(alt_length)
                                else:
                                    vtypes[j].append("INDEL")
                                    vsizes[j].append(len_diff)
                        else:
                            alt_length = 0
                            for idx in alt_idx:
                                alt_length += len(steps[idx])

                            ref_length = 0
                            inv_length = 0
                            isinv = False
                            for idx in ref_idx:
                                if steps[idx][0] != "<":
                                    ref_length += len(reffile.fetch(*steps[idx][1:]))
                                else:
                                    isinv = True
                                    inv_length += len(reffile.fetch(*steps[idx][1:]))

                            if isinv:
                                vtypes[j].append("INV")
                                vsizes[j].append(inv_length)
                                if len(refseq) == ref_length + inv_length + alt_length:
                                    if alt_length > 0:
                                        vtypes[j].append("SNP")
                                        vsizes[j].append(alt_length)
                                else:
                                    if alt_length > 0:
                                        vtypes[j].append("INDEL")
                                        len_diff = (ref_length + inv_length + alt_length) - len(refseq)
                                        vsizes[j].append(len_diff)
                            else:
                                vtypes[j].append("UNKNOWN")
                                vsizes[j].append(-1)


                    else:
                        vtypes[j].append("INDEL")
                        vsizes[j].append(-len(refseq))

                vtype_str = ",".join("|".join(vtype) for vtype in vtypes)
                vsize_str = ",".join("|".join(map(str, vsize)) for vsize in vsizes)
                outfile.write(line)
                outfile.write(f"\t{vtype_str}\t{vsize_str}\n")
