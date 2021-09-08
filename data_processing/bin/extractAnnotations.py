#!/usr/bin/env python3

# Michael VanInsberghe
# Prepare annotations dataframe
# Returns all transcripts in GTF
# Identifies a set of canonical genes

import os, sys, argparse
import pandas as pd
import feather
import gffutils
import pickle

def gene_model(db, tx):
    """ returns the gene model coordinates for a transcript """
    # start codon is first 3 bases of CDS
    # stop codon is last 3 bases of utr3
    model = {'utr5': 0, 'cds': 0, 'utr3': 0}

    if tx.strand == '+':
        childOrder = 'end'
    else:
        childOrder = '-start'

    if 'transcript_type' in tx.attributes:
        if tx.attributes['transcript_type'][0] == 'protein_coding':
            T = db.children(tx, order_by=childOrder, featuretype=['UTR','CDS','start_codon','stop_codon'])
        else:
            T = db.children(tx, order_by=childOrder, featuretype=['exon'])
    elif 'transcript_biotype' in tx.attributes:
        if tx.attributes['transcript_biotype'][0] == 'protein_coding':
            T = db.children(tx, order_by=childOrder, featuretype=['UTR','CDS','start_codon','stop_codon'])
        else:
            T = db.children(tx, order_by=childOrder, featuretype=['exon'])

    current = 'utr5'
    seenCDS = False
    start = []
    end = []
    for child in T:
        #print(child)
        #print(child.featuretype, end="\t")
        if child.featuretype == 'start_codon' or child.featuretype == 'CDS':
            current = 'cds'
            seenCDS = True
        if child.featuretype == 'start_codon':
            continue

        if child.featuretype == 'stop_codon' or (child.featuretype == 'UTR' and seenCDS):
            current = 'utr3'
        if child.featuretype == 'stop_codon':
            continue

        model[current] += child.end - child.start + 1
        start.append(child.start)
        end.append(child.end)

    return {'model': model,
            'start': start,
            'end': end }

parser = argparse.ArgumentParser()
parser.add_argument("--gtf", "-g", help = "path to gtf file")
parser.add_argument("--out", "-o", help = "path to output feather file")
parser.add_argument("--feather", help = "Output reads as feather file") 

args = parser.parse_args()

if args.gtf is None: 
    sys.exit("please provide gtf file")

gtfBase = os.path.splitext( os.path.basename(args.gtf) )[0]
dbOutPath = os.path.join( os.path.dirname(args.gtf), gtfBase + ".annotations.db" )
pickleOutPath = os.path.join( os.path.dirname(args.gtf), gtfBase + ".annotations.pickle" )

if(args.feather is not None):
    featherOutPath = os.path.join( os.path.dirname(args.gtf), gtfBase + ".annotations.feather" )
else:
    csvOutPath = os.path.join( os.path.dirname(args.gtf), gtfBase + ".annotations.csv.gz" )




G = gffutils.create_db(
    args.gtf,
    dbfn=dbOutPath,
    force=True,
    verbose=True,
    merge_strategy='merge',
    disable_infer_transcripts=True,
    disable_infer_genes=True
)

# G = gffutils.FeatureDB(dbOutPath)

annot = []
canonical = {}  # indexed by ENSG, points to annot row
i = 0

for l, feature in enumerate(G.features_of_type('transcript')):
    gene_id = feature.attributes['gene_id'][0]
    
    if 'gene_type' in feature.attributes:
        gene_type = feature.attributes['gene_type'][0]
    elif 'gene_biotype' in feature.attributes:
        gene_type = feature.attributes['gene_biotype'][0]
    else:
        gene_type = "undefined"
    
    if 'transcript_type' in feature.attributes:
        transcript_type = feature.attributes['transcript_type'][0]
    elif 'transcript_biotype' in feature.attributes:
        transcript_type = feature.attributes['transcript_biotype'][0]
    else:
        transcript_type = "undefined"


    gene_info = gene_model(G,feature)
    model = gene_info['model']

    if gene_type == 'protein_coding':
        if 'tag' in feature.attributes:
            # only select protein-coding genes in appris_principal_1 set
            if any('appris_principal' in tag for tag in feature.attributes['tag']):
                if gene_id in canonical:
                    # if exists, replace if the CDS is longer
                    # or if CDS are the same, if the transcripe length is longer
                    if model['cds'] > annot[ canonical[gene_id] ]['l_cds']:
                        canonical[gene_id] = l
                    elif model['cds'] == annot[ canonical[gene_id] ]['l_cds'] and sum(model.values()) > annot[ canonical[gene_id] ]['l_tr']:
                        canonical[gene_id] = l
                else:
                    canonical[gene_id] = l 
    else:
        if gene_id in canonical:
            # if exists, replace if transcript is longer
            if sum(model.values()) > annot[ canonical[gene_id] ]['l_tr']:
                canonical[gene_id] = l
            else:
                canonical[gene_id] = l
    annot.append(
            {   'transcript_id': feature.attributes['transcript_id'][0],
                'gene_id': feature.attributes['gene_id'][0],
                'gene_name': feature.attributes['gene_name'][0],
                'gene_type': gene_type, 
                'transcript_type': transcript_type,
                'chr': feature.seqid,
                'strand': feature.strand,
                'l_tr': sum(model.values()),
                'l_utr5': model['utr5'],
                'l_cds': model['cds'],
                'l_utr3': model['utr3'],
                'set': 'none',
                'start': gene_info['start'],
                'end': gene_info['end']
            })
    #print(annot[i])
    i += 1

for g in canonical:
    annot[ canonical[g] ]['set'] = "canonical"

# remove start and end coordinates for R export
annotations = [{k: v for k, v in d.items() if k != 'start' if k != 'end' } for d in annot]

if(args.feather is not None):
    feather.write_dataframe(pd.DataFrame(annotations), featherOutPath)
else:
    pd.DataFrame(annotations).to_csv(csvOutPath, mode = "w", compression = "infer", index = False)

annodict = {}
for item in annot:
    name = item.pop('transcript_id')
    annodict[name] = item

#print(annodict)
with open(pickleOutPath, 'wb') as f:
    pickle.dump(annodict, f, pickle.HIGHEST_PROTOCOL)
