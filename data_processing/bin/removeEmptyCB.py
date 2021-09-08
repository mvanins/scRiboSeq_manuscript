#!/usr/bin/env python3
# Michael VanInsberghe 2021-02-25
# removes CB:Z:- and UB:Z:- tags added by STAR

import os, sys, argparse
import pysam
import csv
import math

argparser = argparse.ArgumentParser()
argparser.add_argument("--bam", "-b", help = "Input bam file")
argparser.add_argument("--out", "-o", help = "Output bam file. Default appends _CB to filename. E.g., Input_CB.bam")
argparser.add_argument("--threads", "-t", help = "Number of threads", type=int)

args = argparser.parse_args()

if args.bam is None:
    sys.exit("please provide bam file")

if args.out is None:
    bamName = os.path.basename(args.bam) # bam filename
    bamBase = os.path.splitext(bamName)[0] # remove extension
    outPath = os.path.join ( os.path.dirname(args.bam), bamBase + "_stripped.bam" )
else:
    outPath = args.out

if args.threads is None:
    pysthr = 1
else:
    pysthr = max( 1, math.floor( (args.threads - 1)/2 ) )


# process bam
try:
    bam = pysam.AlignmentFile(args.bam, "rb", threads=pysthr)
except OSError:
    print("Could not open input file:", args.bam)
    sys.exit()

try:
    outbam = pysam.AlignmentFile(outPath, "wb", template = bam, threads=pysthr)
except OSError:
    print("Could not open output file:", outPath)


for read in bam:
    try:
        CB = read.get_tag("CB")
    except KeyError:
        continue

    try:
        UB = read.get_tag("UB")
    except KeyError:
        continue

    if str(CB) == '-':
        read.set_tag('CB', None)

    if str(UB) == '-':
        read.set_tag('UB', None)

    outbam.write(read)

bam.close()
outbam.close()

