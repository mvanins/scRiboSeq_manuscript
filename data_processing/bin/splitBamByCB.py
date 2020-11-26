#!/usr/bin/env python3
# Michael VanInsberghe 2020-04-02
# splits BAM files based on supplied tag CB default


import os, sys
import pysam
import csv
import argparse
import math

argparser = argparse.ArgumentParser()
argparser.add_argument("--bam","-b", help="set bam file")
argparser.add_argument("--CB","-c", help="tsv of barcodes in barcode id well format")
argparser.add_argument("--flag", "-f", help="BAM tag containing cell barcode", default="CB")
argparser.add_argument("--threads", "-t", help="threads for decompressing input BAM file", type=int)

args = argparser.parse_args()

if args.bam is None or args.CB is None:
    sys.exit("please provide both bam and CB")

if args.threads is None:
    pysthr = 1
else:
    pysthr = max( 1, math.floor( (args.threads - 1)/2 ) )
    #pysthr = max( 1, args.threads-1)


# barcodes = {}
# with open(args.CB) as barcodeFile:
#     csvReader = csv.reader(barcodeFile, delimiter='\t')
#     for row in csvReader:
#         barcodes[row[0]] = row[1]
# 
# barcodeFile.close()

barcodes = {}
with open(args.CB) as barcodeFile:
    for row in barcodeFile:
        row = row.rstrip("\r\n")
        if not row.strip():
            continue

        bcd = row.split("\t")[0]
        #barcodes[bcd] = row.split("\t")[1]
        barcodes[bcd] = 1



try:
    bam = pysam.AlignmentFile(args.bam, "rb", threads = pysthr)
except:
    sys.exit("bam file %s not found" % args.bam)


# put output files in subfolder based on library name
bamname = os.path.basename(args.bam) # filename
bambase = os.path.splitext(bamname)[0] # remove extension

if "_" in bambase:
    library = bambase.split("_")[0]
    extras = "_" + bambase.split("_",1)[1]
else:
    library = bambase
    extras = ""

dirname = os.path.dirname(args.bam)
outdir = os.path.join( os.path.dirname(args.bam) , bambase )

if not os.path.exists(outdir):
    os.mkdir(outdir)

# open output bam files
outbams = {}
for barcode in barcodes:
    outbams[barcode] = pysam.AlignmentFile( os.path.join( outdir, library + "_" + barcode + extras + ".bam" ), "wb", template = bam , threads = pysthr)


# iterate through bam file, segregate reads based on given tag
nocb = 0
nreads = 0
for read in bam:
    nreads += 1
    try:
        CB = read.get_tag(args.flag)
    except KeyError:
        nocb += 1
        continue

    if CB in barcodes:
        outbams[CB].write(read)
    else:
        continue

bam.close()
for bamout in outbams.values():
    bamout.close()
print(f"Warning: {args.flag} tag not found in {nocb} of {nreads} reads")
