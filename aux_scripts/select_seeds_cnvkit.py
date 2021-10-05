#!/usr/bin/env python
import sys
import os
import threading
from subprocess import call
import argparse
import gzip

## code adapted from  PrepareAA.py commit 38e758c19776f02cd0e920bde2513536cccc489f 
## https://github.com/jluebeck/PrepareAA


# Read the CNVkit .cns files
def convert_cnvkit_cnv_to_seeds(cnvkit_output_directory, bam):
    base = os.path.splitext(os.path.basename(bam))[0]
    with open(cnvkit_output_directory + base + ".segment.cns") as infile, open(cnvkit_output_directory + base + "_CNV_GAIN.bed",
                                                                       'w') as outfile:
        head = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            s, e = int(fields[1]), int(fields[2])
            cn_r = float(fields[4])
            cn = 2 ** (cn_r + 1)
            if cn >= args.cngain and e - s >= args.cnsize_min:
                outline = "\t".join(fields[0:3] + ["CNVkit", str(cn)]) + "\n"
                outfile.write(outline)

    return cnvkit_output_directory + base + "_CNV_GAIN.bed"

# MAIN #
if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(
        description="A simple pipeline wrapper for AmpliconArchitect, invoking alignment, variant calling, "
                    "and CNV calling prior to AA. The CNV calling is necesary for running AA")
    #parser.add_argument("-s", "--sample_name", help="sample name", required=True)
    parser.add_argument("-o", "--output_directory", help="output directory names (will create if not already created)", required=True)
    parser.add_argument("--sorted_bam", help="Sorted BAM file (aligned to an AA-supported reference.)", required=True)
    parser.add_argument("--cngain", type=float, help="CN gain threshold to consider for AA seeding", default=4.5)
    parser.add_argument("--cnsize_min", type=int, help="CN interval size (in bp) to consider for AA seeding",
                        default=50000)
    args = parser.parse_args()
    args.cnv_bed = convert_cnvkit_cnv_to_seeds(args.output_directory + "/", args.sorted_bam)
