#!/usr/bin/env python3

import argparse, logging, tarfile, subprocess, os, sys
from datetime import datetime
from collections import Counter
import pandas as pd
from Bio import SeqIO
import tempfile

parser = argparse.ArgumentParser(description="Adds taxon lineage to DIAMOND BlastX output")

parser.add_argument('-i',
                    '--infile',
                    help="Input tsv file",
                    type=str,
                    required = True)

parser.add_argument('-c',
                    '--contig_file',
                    help="Input contig file",
                    type=str,
                    required = True)

parser.add_argument('-o',
                    '--outfile',
                    help="Output tsv file",
                    type=str,
                    required = True)

parser.add_argument('-u',
                    '--unannotated',
                    help="Output tsv file",
                    type=str,
                    required = True)

parser.add_argument('-log',
                    '--logfile',
                    help="Name of logfile",
                    default = datetime.now().strftime("logfile_%d-%m-%Y.log"),
                    type=str,
                    required = False)

parser.add_argument('--quiet',
                    action='store_true',
                    default = False,
                    required = False)

def create_taxon_dict():
    '''Download and create taxon dictionary'''
    logging.info("Downloading new_taxdump..")
  
    _, temp = tempfile.mkstemp()

    subprocess.run(['rsync','rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz',temp], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    
    logging.info("Creating taxon dictionary..")

    taxon_dict = {}
    with tarfile.open(temp, "r:gz") as tar:
        for line in tar.extractfile("rankedlineage.dmp"):
            line = line.decode().strip().replace('\t','').split('|')
            taxon_dict[line[0]] = line[1:-1]

    #Add merged entries
    with tarfile.open(temp, "r:gz") as tar:
        for line in tar.extractfile("merged.dmp"):
            line = line.decode().strip().replace('\t','').split('|')
            taxon_dict[line[0]] = taxon_dict[line[1]]

    #Add deleted entries
    with tarfile.open(temp, "r:gz") as tar:
        for line in tar.extractfile("delnodes.dmp"):
            line = line.decode().strip().replace('\t','').split('|')
            taxon_dict[line[0]] = taxon_dict['1']
    
    os.unlink(temp)
    
    return(taxon_dict)

class contig_annotation:

    def __init__(self, contig_id):
        self.contig_id = contig_id
        self.annotation_list = []
        
    def get_lineage(self, tax_id):
        try:
            return(taxon_dict[tax_id])
        except Exception as e:
            return(["unknown"]+["" for x in range(9)])

    def add_annotation(self, line):

        if len(line) == 13:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, staxids = line
            staxids = staxids.split(';')
        elif len(line) == 12:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line
            staxids = [1]
        else:
            return

        if len(staxids) == 1:
            staxids = staxids[0]
        else:
            majority_taxid = sorted(Counter(staxids))[-1]
            staxids = majority_taxid

        lineage = self.get_lineage(staxids)

        if len(lineage) < 9:
            print(f"DEBUG: Short lineage for TaxID '{staxids}' (Contig: {self.contig_id})."
                  f"Lineage ({len(lineage)} ranks): {lineage}", file = sys.stderr)

        lineage.extend([''] * (9 - len(lineage)))

        new_annotation = {
            'qseqid': qseqid,
            'sseqid': sseqid,
            'pident': pident,
            'length': length,
            'mismatch': mismatch,
            'gapopen': gapopen,
            'qstart': qstart,
            'qend': qend,
            'sstart': sstart,
            'send': send,
            'evalue': evalue,
            'bitscore': bitscore,
            'staxids': staxids,
            'taxon_name': lineage[0],
            'species': lineage[1],
            'genus': lineage[2],
            'family': lineage[3],
            'order': lineage[4],
            'class': lineage[5],
            'phylum': lineage[6],
            'kingdom': lineage[7],
            'superkingdom': lineage[8]
        }

        self.annotation_list.append(new_annotation)

    def get_best_hit(self):
        if len(self.annotation_list) == 0:
            return(None)
        annotation_df = pd.DataFrame(self.annotation_list)

        annotation_df['order_column'] = annotation_df['bitscore'].astype(float) / annotation_df['length'].astype(float)
        best_hit = annotation_df.sort_values(by='order_column', ascending=False).iloc[0].drop(labels=['order_column']).astype(str)
        return(best_hit.tolist())

    def format(self):
        best_hit = self.get_best_hit()
        if best_hit == None:
            return(None)
        return('\t'.join(best_hit))
        
def parse_diamond(infile, outfile, total):
    logging.info("Creating annotated file")
    annotated_set = set()
    progress = 1
    with open(infile,"r") as infile, open(outfile,"w") as outfile:
        contig = contig_annotation(None)
        for line in infile:
            line = line.strip()
            line_split = line.split('\t')

            contig_id = line_split[0]

            annotated_set.add(contig_id)
            
            if contig_id == contig.contig_id:
                contig.add_annotation(line_split)
                continue
            else:
                progress += 1
                if not args.quiet:
                    print(f'Progress: {progress}/{total}',end="\r")
                if contig.contig_id != None:
                    output = contig.format()
                    if output != None:
                        print(contig.format(), file = outfile)
                contig = contig_annotation(contig_id)
                contig.add_annotation(line_split)
    return(annotated_set)

if __name__ == "__main__":
    args = parser.parse_args()
        
    if args.quiet:
        logging.basicConfig(handlers=[
                                logging.FileHandler(args.logfile),
                                logging.StreamHandler()
                            ],
                            level=logging.CRITICAL,
                            format='%(levelname)s %(asctime)s %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(handlers=[
                                logging.FileHandler(args.logfile),
                                logging.StreamHandler()
                            ],
                            level=logging.INFO,
                            format='%(levelname)s %(asctime)s %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')

    global taxon_dict

    taxon_dict = create_taxon_dict()

    logging.info("Creating contig set")
    contig_set = set()
    for record in SeqIO.parse(args.contig_file, "fasta"):
        contig_set.add(record.id)

    total = len(contig_set)

    annotated_set = parse_diamond(args.infile, args.outfile, total)
    
    logging.info("Creating unannotated file")
    unannotated = contig_set.difference(annotated_set)

    with open(args.unannotated, "w") as f:
        for contig in unannotated:
            print(contig, file = f)
