#!/usr/bin/env python3
# Michael VanInsberghe 2020-04-08
# Parse bam file, extracting read information 

import os, sys, argparse
import pickle, feather
import pandas as pd
import numpy as np
import gffutils
import pysam
import datetime
import uuid

def reverse_complement(seq):
    """ returns the reverse complement of a sequence """

    trans = seq.maketrans("ABCDGHMNRSTUVWXYabcdghmnrstuvwxy", "TVGHCDKNYSAABWXRtvghcdknysaabwxr")
    return seq.translate(trans)[::-1]

def txToGenomeList_(transInfo, coords = []):
    """ returns the 0-based genome coordinates for a given trancript_id and a list of coordinates """
    
    if transInfo['strand'] == '+':
        childOrder = 'end'
    else:
        childOrder = '-start'
    
    coord_rel = coords.copy()
    genomic = [None]*len(coords)
    
    for start, end in zip(transInfo['start'], transInfo['end']):
        exon_length = end - start + 1
        
        hits = [i for i, l in enumerate(coord_rel) if l-exon_length <= 0]
        if hits:
            # at least one hit is in this exon
            if transInfo['strand']  == '+':
                replacements = [start + coord_rel[i] for i in hits]
            else:
                replacements = [end - coord_rel[i] for i in hits]
            
            for (index, replacement) in zip(hits,replacements):
                genomic[index] = replacement
                coord_rel[index] = float('inf') 

            if all(v is not None for v in genomic): 
                return {'chr': transInfo['chr'], 'genome': genomic }	
    
        coord_rel[:] = [x - exon_length for x in coord_rel]
    
    return {'chr': transInfo['chr'], 'genome': genomic }



def fetchBaseList_(transInfo,fasta,coords = []):
    """ fetches bases for given transcript coordinates """

    genomic = txToGenomeList_(transInfo, coords)

    if genomic is not None:
        bases = [fasta[genomic['chr']][x-1] if (x is not None and x <= len(fasta[genomic['chr']])-1)  else "" for x in genomic['genome']]
        #bases = [fasta.fetch(genomic['chr'], x-1, x) if x is not None else "" for x in genomic['genome'] ]
        if transInfo['strand'] == '-':
            bases[:] = [reverse_complement(base) for base in bases]

        return bases

    else:
        return [""]*len(coords)

def cut_bases_(transInfo,fasta,read):
    """ returns the bases around the 5' and 3' end of a read """

    start = read.reference_start 
    end = read.reference_end - 1 

    bases = fetchBaseList_(transInfo,fasta,[start, start-1, start-2, start-3, start-4, start-5, start-6, start-7, start-8,
        start+1, start+2, start+3, start+4, start+5, start+6, start+7,start+8,
        end,end-1,end-2,end-3,end-4,end-5,end-6,end-7,end-8,
        end+1,end+2,end+3,end+4,end+5,end+6,end+7,end+8])
    
    return(dict(zip(['base5','base5m1','base5m2','base5m3','base5m4','base5m5','base5m6','base5m7','base5m8',
        'base5p1','base5p2','base5p3','base5p4','base5p5','base5p6','base5p7','base5p8',
        'base3','base3m1','base3m2','base3m3','base3m4','base3m5','base3m6','base3m7','base3m8',
        'base3p1','base3p2','base3p3','base3p4','base3p5','base3p6','base3p7','base3p8'],
        bases)))

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




parser = argparse.ArgumentParser()
parser.add_argument("--bam", "-b", help = "path to input bam file")
parser.add_argument("--trans", "-t", help = "path to input transcriptome.pickle model annotations from extractAnnotations")
parser.add_argument("--out", "-o", help = "path to output feather file for read information")
parser.add_argument("--fasta", "-f", help = "path to input fasta file. fai required")
parser.add_argument("--CB", "-c", help = "BAM tag containing cell barcode", default = "CB")
parser.add_argument("--set", "-s", help = "Only include transcripts in supplied set")
parser.add_argument("--feather", help = "Output reads as feather file. Requires substantially more memory", action='store_const', const=1) 
parser.add_argument("--id", help = "Output read ID. Cannot be combined with --feather", action='store_const', const=1)

args = parser.parse_args()

if args.trans is None or args.bam is None or args.fasta is None:
    sys.exit("Please provide: --(t)rans --(b)am --(f)asta")

if args.feather is not None and args.id is not None:
    sys.exit("--id and --feather cannot be combined")

# output file
if args.out is None:
    bamBase = os.path.splitext( os.path.basename(args.bam) )[0]
    if(args.feather is not None):
        outPath = os.path.join(bamBase + ".feather")
    else:
        outPath = os.path.join(bamBase + ".csv.gz")
else:
    outPath = args.out

if(args.feather is not None):
    try:
        tempfn = str(uuid.uuid4()) + ".pkl"
        tempfh = open(tempfn,'ab')
    except:
        print(f"could not open temporary output file {tempfn}")
        sys.exit()
#else:
#    try:
#        csvfh = open(outPath, 'wb') 
#    except:
#        print(f"Could not open output file {outPath}")
#        sys.exit()


# annotations
with open(args.trans, 'rb') as f:
    print("loading annotations")
    annot = pickle.load(f)

if(args.set is not None):
    print(f"Only keeping annotations in {args.set} set")
    nk = len(annot.keys())
    print(f"length of input annotations {nk}")
        
    for k in list(annot.keys()):
        if annot[k]['set'] != args.set:
            del annot[k]
        
    nk = len(annot.keys())
    print(f"length of reduced annotations {nk}")

# genome
fasta = loadGenome(args.fasta)

# input bam file
try:
    bam = pysam.AlignmentFile(args.bam, "rb", threads = 2)
except OSError:
    print("could not open input file:", args.bam)
    sys.exit()

CBcount = 0
readcount = 0
for read in bam.head(100):
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

print(f"Collect cell barcodes from BAM tag {args.CB}: {CollectCB}")
print(datetime.datetime.now(), end = "\t")
print(f"processsing {args.bam}")
i = 0
ilast = 0
tempreads = []
c = []
for read in bam.fetch():

    if (i>ilast) & (i % 200000 == 0):
        print(datetime.datetime.now(), end = "\t")
        print(f"processed {i} alignments")
        
        reads = pd.DataFrame(tempreads)
        if(args.feather is not None):
            pickle.dump(reads,tempfh)
            c.append(len(reads))
        else:
            if(ilast == 0):
                reads.to_csv(outPath, mode = "w", compression = "infer", index = False)
            else:
                reads.to_csv(outPath, header = None, mode = "a", compression = "infer", index = False)

        tempreads.clear()
        ilast = i

    if read.reference_name in annot:
        transcript_id = read.reference_name
	
        bases = cut_bases_(annot[transcript_id], fasta, read)


        if annot[transcript_id]['transcript_type'] == 'protein_coding':
            cds_start = annot[transcript_id]['l_utr5']
            cds_end = annot[transcript_id]['l_utr5'] + annot[transcript_id]['l_cds']
        else:
            cds_start = 0
            cds_end = annot[transcript_id]['l_tr']

        strand = "+"
        if read.is_reverse:
            strand = "-"

        readinfo = {
                'transcript_id': read.reference_name,
                'cut5': read.reference_start,
                'cut3': read.reference_end,
                'length': read.query_alignment_length,
                'read_length': read.query_length,
                'read_strand': strand,
                'cds_start': cds_start,
                'cds_end': cds_end,
                'frame5': (read.reference_start - cds_start) % 3,
                'frame3': (read.reference_end - cds_start -1 ) % 3
                }

        if args.id is not None:
            readinfo['id'] = read.query_name

        if CollectCB:
            try:
                readinfo[args.CB] = read.get_tag(args.CB)
            except KeyError:
                readinfo[args.CB] = ""

        tempreads.append({**readinfo, **bases})

    else:
        continue

    #print(reads[i])
    i += 1

bam.close()

print(datetime.datetime.now(), end = "\t")
print(f"End of bam file. processed {i} alignments")

reads = pd.DataFrame(tempreads)
if(args.feather is not None):
    pickle.dump(reads,tempfh)
    c.append(len(reads))
    ncols = len(tempreads[0])
    colnames = reads.columns
    tempfh.close()
else:
    if(ilast == 0):
        reads.to_csv(outPath, mode = "w", compression = "infer", index = False)
    else:
        reads.to_csv(outPath, header = None, mode = "a", compression = "infer", index = False)
    #csvfh.close()

del tempreads, reads, fasta


print("Finished!")

if(args.feather is not None):
    print(f"Merging and writing to output file {outPath}")
    
    with open(tempfn, 'rb') as tempfh:
        reads_all = pickle.load(tempfh)
        offset = len(reads_all)
        reads_all = reads_all.append(pd.DataFrame(np.empty(sum(c[1:])*ncols).reshape(-1,ncols), columns = colnames))
        #reads_all = reads_all.append(pd.DataFrame(np.empty(sum(c[1:])*4).reshape(-1,4)))
    
        for size in c[1:]:
            df=pickle.load(tempfh)
            reads_all.iloc[offset:offset+size]=df.values 
            offset+=size
    
    os.remove(tempfn)
    feather.write_dataframe(reads_all, outPath)
