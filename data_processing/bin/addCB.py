#!/usr/bin/env python3
# Michael VanInsberghe 2020-04-02
# correct raw cell barcode CR to CB with a single-base hamming expansion
# add CB to bam file
# divides given threads between (de)compression threads for pysam read/write


import os, sys, argparse
import pysam
import csv
import math

argparser = argparse.ArgumentParser()
argparser.add_argument("--bam", "-b", help = "Input bam file")
argparser.add_argument("--CB", "-c", help = "tsv of barcodes with barcode in first column")
argparser.add_argument("--out", "-o", help = "Output bam file. Default appends _CB to filename. E.g., Input_CB.bam")
argparser.add_argument("--threads", "-t", help = "Number of threads", type=int)

args = argparser.parse_args()

if args.bam is None or args.CB is None:
    sys.exit("please provide both bam and CB")

if args.out is None:
    bamName = os.path.basename(args.bam) # bam filename
    bamBase = os.path.splitext(bamName)[0] # remove extension
    outPath = os.path.join ( os.path.dirname(args.bam), bamBase + "_CB.bam" )
else:
    outPath = args.out

if args.threads is None:
    pysthr = 1
else:
    pysthr = max( 1, math.floor( (args.threads - 1)/2 ) )


# read in barcodes and expand 1 hamming distance
bases = "ATGCN"
corrbcd = {}

with open(args.CB, "r") as fh:
    for row in fh:
        row = row.rstrip("\r\n")
        if not row.strip():
            continue	
        
        bcd = row.split("\t")[0]
        corrbcd[bcd] = bcd
        for p in range( len(bcd) ):
            alternatives = bases.replace(bcd[p],"")
            for m in alternatives:
                mutated = list(bcd)
                mutated[p] = m;
                corrbcd["".join(mutated)] = bcd


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
        CR = read.get_tag("CR")
    except KeyError:
        continue

    if CR in corrbcd:
        read.set_tag('CB',corrbcd[CR])

    outbam.write(read)

bam.close()
outbam.close()

