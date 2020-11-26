#!/usr/bin/env python3
# Michael VanInsberghe 2020-04-30
# Extracts codon table for protein-coding transcripts 

import os, sys, argparse, warnings
import pickle, feather
import pandas as pd
import numpy as np
import pysam
# import gffutils
# import datetime
# import uuid

translation = {'TTT': {'single': 'F', 'three': 'Phe'},
        'TTC': {'single': 'F', 'three': 'Phe'},
        'TTA': {'single': 'L', 'three': 'Leu'},
        'TTG': {'single': 'L', 'three': 'Leu'},
        'CTT': {'single': 'L', 'three': 'Leu'},
        'CTC': {'single': 'L', 'three': 'Leu'},
        'CTA': {'single': 'L', 'three': 'Leu'},
        'CTG': {'single': 'L', 'three': 'Leu'},
        'ATT': {'single': 'I', 'three': 'Ile'},
        'ATC': {'single': 'I', 'three': 'Ile'},
        'ATA': {'single': 'I', 'three': 'Ile'},
        'ATG': {'single': 'M', 'three': 'Met'},
        'GTT': {'single': 'V', 'three': 'Val'},
        'GTC': {'single': 'V', 'three': 'Val'},
        'GTA': {'single': 'V', 'three': 'Val'},
        'GTG': {'single': 'V', 'three': 'Val'},
        'TCT': {'single': 'S', 'three': 'Ser'},
        'TCC': {'single': 'S', 'three': 'Ser'},
        'TCA': {'single': 'S', 'three': 'Ser'},
        'TCG': {'single': 'S', 'three': 'Ser'},
        'CCT': {'single': 'P', 'three': 'Pro'},
        'CCC': {'single': 'P', 'three': 'Pro'},
        'CCA': {'single': 'P', 'three': 'Pro'},
        'CCG': {'single': 'P', 'three': 'Pro'},
        'ACT': {'single': 'T', 'three': 'Thr'},
        'ACC': {'single': 'T', 'three': 'Thr'},
        'ACA': {'single': 'T', 'three': 'Thr'},
        'ACG': {'single': 'T', 'three': 'Thr'},
        'GCT': {'single': 'A', 'three': 'Ala'},
        'GCC': {'single': 'A', 'three': 'Ala'},
        'GCA': {'single': 'A', 'three': 'Ala'},
        'GCG': {'single': 'A', 'three': 'Ala'},
        'TAT': {'single': 'Y', 'three': 'Tyr'},
        'TAC': {'single': 'Y', 'three': 'Tyr'},
        'TAA': {'single': '*', 'three': 'Ter'},
        'TAG': {'single': '*', 'three': 'Ter'},
        'CAT': {'single': 'H', 'three': 'His'},
        'CAC': {'single': 'H', 'three': 'His'},
        'CAA': {'single': 'Q', 'three': 'Gln'},
        'CAG': {'single': 'Q', 'three': 'Gln'},
        'AAT': {'single': 'N', 'three': 'Asn'},
        'AAC': {'single': 'N', 'three': 'Asn'},
        'AAA': {'single': 'K', 'three': 'Lys'},
        'AAG': {'single': 'K', 'three': 'Lys'},
        'GAT': {'single': 'D', 'three': 'Asp'},
        'GAC': {'single': 'D', 'three': 'Asp'},
        'GAA': {'single': 'E', 'three': 'Glu'},
        'GAG': {'single': 'E', 'three': 'Glu'},
        'TGT': {'single': 'C', 'three': 'Cys'},
        'TGC': {'single': 'C', 'three': 'Cys'},
        'TGA': {'single': '*', 'three': 'Ter'},
        'TGG': {'single': 'W', 'three': 'Trp'},
        'CGT': {'single': 'R', 'three': 'Arg'},
        'CGC': {'single': 'R', 'three': 'Arg'},
        'CGA': {'single': 'R', 'three': 'Arg'},
        'CGG': {'single': 'R', 'three': 'Arg'},
        'AGT': {'single': 'S', 'three': 'Ser'},
        'AGC': {'single': 'S', 'three': 'Ser'},
        'AGA': {'single': 'R', 'three': 'Arg'},
        'AGG': {'single': 'R', 'three': 'Arg'},
        'GGT': {'single': 'G', 'three': 'Gly'},
        'GGC': {'single': 'G', 'three': 'Gly'},
        'GGA': {'single': 'G', 'three': 'Gly'},
        'GGG': {'single': 'G', 'three': 'Gly'}}




def reverse_complement(seq):
    """ returns the reverse complement of a sequence """
    
    trans = seq.maketrans("ABCDGHMNRSTUVWXYabcdghmnrstuvwxy", "TVGHCDKNYSAABWXRtvghcdknysaabwxr")
    return seq.translate(trans)[::-1]

def loadGenome(fastafn):
    """ loads the genome into a dict """

    try:
        fasta = pysam.FastaFile(fastafn)
    except IOError:
        print("Could not open input genome file:", fastafn)
        sys.exit()

    print("loading genome sequences")

    sequences = {}
    for name in fasta.references:
        sequences[name] = fasta.fetch(name)

    fasta.close()
    return sequences


def transcriptSeq(transInfo, fasta):
    """ assembles the full exonic transcript sequence given a transcript model """

    seq = ""

    for s, e in zip(transInfo['start'], transInfo['end']):
        exon = fasta[transInfo['chr']][s-1:e:1]
        if(transInfo['strand'] == "-"):
            exon = reverse_complement(exon)
        #print(f"{s}\t{e}\t{exon}")

        seq += exon

    return seq

def codingSeq(transInfo, transeq):
    """ subsets the coding sequence from a transcript sequence """
    coding = transeq[transInfo['l_utr5']:(transInfo['l_utr5'] + transInfo['l_cds'])]
    return coding

def translate(codon,translation):
    """ returns a translation for a codon """

    if codon in translation:
        return translation[codon]

    if "N" in codon:
        # check if it's degenerate
        codons = [codon.replace("N",x) for x in ["A","T","G","C"]]
        AAs = [translation[c]['single'] for c in codons]
        if AAs.count(AAs[0]) == 4:
            return translation[codons[0]]

    return {'single': '', 'three': ''}



parser = argparse.ArgumentParser()
parser.add_argument("--trans", "-t", help = "path to input transcriptome.pickle model annotations from extractAnnotations")
parser.add_argument("--out", "-o", help = "path to output feather file for read information")
parser.add_argument("--fasta", "-f", help = "path to input genome fasta file. fai required")
parser.add_argument("--feather", help = "Output reads as feather file") 


args = parser.parse_args()

if args.trans is None or args.fasta is None:
    sys.exit("Please provide: --(t)rans --(f)asta")

if args.out is None:
    transBase = os.path.splitext( os.path.basename(args.trans) )[0]
    if(args.feather is not None):
        outPath = os.path.join(transBase + ".transcript_codons.feather")
    else:
        outPath = os.path.join(transBase + ".transcript_codons.csv.gz")
else:
    outPath = args.out

# load annotations
with open(args.trans, 'rb') as f:
    print("loading annotations")
    annot = pickle.load(f)

# load genome
fasta = loadGenome(args.fasta)

codons = []
for t in annot:
    transInfo = annot[t]
    
    if transInfo['transcript_type'] != "protein_coding":
        continue
    
    coding = codingSeq(transInfo, transcriptSeq(transInfo, fasta))
#    if transInfo['set'] == "canonical" and (len(coding) % 3) != 0:
#        warnings.warn(f"Warning: {t} in canonical set skipped because length not multple of 3")
#        continue
#    elif (len(coding) % 3) != 0:
#        continue

    start = transInfo['l_utr5']
    #print(t, end="\t")
    #print(len(coding) % 3)
    for codon in [coding[i:i+3] for i in range(0, len(coding), 3)]:
        codons.append({'transcript_id': t,
            'codon': codon,
            'single': translate(codon,translation)['single'],
            'three': translate(codon,translation)['three'],
            'codon_position': start})
        start = start + 3

if(args.feather is not None):
    feather.write_dataframe(pd.DataFrame(codons), outPath)
else:
    pd.DataFrame(codons).to_csv(outPath, mode = "w", compression = "infer", index = False)

