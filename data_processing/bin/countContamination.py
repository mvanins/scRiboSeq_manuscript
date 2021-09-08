#!/usr/bin/env python3
# Michael VanInsberghe 2021-03-02
# Counts reads mapping to contaminants, segregating either as 
#   - tRNAs: contig beginning with tRNACluster
#   - ribosomal: contig beginning with anything else
# if read maps to both tRNA and rRNA, assumed to stem from tRNA
# inflexible, but good enough for now

import os, sys, argparse
import pysam
import csv
import math
from collections import defaultdict
import pandas as pd

argparser = argparse.ArgumentParser()
argparser.add_argument("--bam", "-b", help = "Input bam file")
argparser.add_argument("--name", "-n", help = "Sample name. Default takes entire name before first underscore")
argparser.add_argument("--CB", "-c", help = "BAM tag containing cell barcode", default = "CB")
argparser.add_argument("--out", "-o", help = "Output count csv file. Default appends _contamination_counts.csv.gz to filename. E.g., Input_contamination.csv.gz")
argparser.add_argument("--threads", "-t", help = "Number of threads", type=int)

args = argparser.parse_args()

if args.bam is None:
    sys.exit("please provide input bam file in --bam")

if args.out is None:
    bamName = os.path.basename(args.bam) # bam filename
    bamBase = os.path.splitext(bamName)[0] # remove extension
    outPath = os.path.join ( os.path.dirname(args.bam), bamBase + "_contamination_counts.csv.gz" )
else:
    outPath = args.out

if args.threads is None:
    pysthr = 1
else:
    pysthr = max( 1, math.floor( (args.threads - 1)/2 ) )

if args.name is None:
    bamName = os.path.basename(args.bam) # bam filename
    name = bamName.split('_',1)[0]
else:
    name = args.name

# read bam
try:
    bam = pysam.AlignmentFile(args.bam, "rb", threads=pysthr)
except OSError:
    print("Could not open input file:", args.bam)
    sys.exit()

# check if bulk
CBcount = 0
readcount = 0
for read in bam.head(200):
    readcount += 1
    try:
        CB = read.get_tag(args.CB)
        CBcount += 1
    except KeyError:
        pass
if CBcount/readcount > .5:
    CollectCB = True
else:
    CollectCB = False


classCounts = defaultdict(lambda: defaultdict(int))

for read in bam:

    if CollectCB:
        # skip reads without assigned CB
        try:
            CB = read.get_tag(args.CB)
        except KeyError:
            continue
    else:
        CB = name

    readid = read.query_name
    contig = read.reference_name

    if contig is None:
        classCounts[CB]['Unaligned'] += 1
        continue

    if contig.startswith('tRNA'):
        classCounts[CB]['tRNA'] += 1
        continue
    else:
        classCounts[CB]['rRNA'] += 1

bam.close()

counts = pd.DataFrame(classCounts).T
counts.to_csv(outPath, mode = "w", compression = "infer", index = True)
