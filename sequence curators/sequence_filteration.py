#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import argparse
import time
from Bio import SeqIO
from Bio.Seq import reverse_complement

def parsing_filter(remove_file):                      #if only names not fasta files
    with open(remove_file, 'r')as f:
        all_lines=f.readlines()
        if '>' in all_lines[1]:
            names_to_del=[line.rstrip()[1:] for line in all_lines if line != '']
        else:
            names_to_del=[line.rstrip() for line in all_lines if line != '']
    return names_to_del

def remover_by_names (start_file,names_to_del):
    records=SeqIO.parse(start_file,'fasta')
    curated=[record for record in records if record.description not in names_to_del]
    SeqIO.write(curated,end_file, "fasta")

def remover_by_seq (start_file,remove_file,database):
    records=SeqIO.parse(start_file,'fasta')
    sequences_to_del=[str(record.seq) for record in SeqIO.parse(remove_file,'fasta')]
    if database == 'n':
        sequences_to_del+=[reverse_complement(seq) for seq in sequences_to_del]
    curated=[record for record in records if record.seq not in sequences_to_del]
    SeqIO.write(curated,end_file, "fasta")

                  #begining of the code !!!

parser = argparse.ArgumentParser(prog='Filtering program',usage='\n%(prog)s :filter nucleotide and/or protein database from a list of names or sequences (by exact match).\n\n Eslam S.Ibrahim\n\n eslam.ebrahim@pharma.cu.edu.eg')
parser.add_argument('-out',dest='end_file',required=True)
parser.add_argument('-in',dest='start_file',type=argparse.FileType('r'),required=True)
parser.add_argument('-filter',dest='remove_file',required=True)
parser.add_argument('-flt_mod',dest='filter_mode',choices=['seq','name'],default='seq',help='if you want to filter out sequences by names only',required=True)
parser.add_argument('-p',dest='database',action='store_const',const='p',help='protein sequences')
parser.add_argument('-n',dest='database',action='store_const',const='n',help='nucleotide sequences')
args=parser.parse_args()

start_file=args.start_file
remove_file=args.remove_file
end_file=args.end_file
filter_mode=args.filter_mode
database=args.database

start_time = time.clock()

if filter_mode=='name':
    remover_by_names(start_file,parsing_filter(remove_file))

elif filter_mode=='seq':
    remover_by_seq(start_file,remove_file,database)


time_of_calc = time.clock() - start_time
print (time_of_calc, "seconds")                                      #in case if we want the time
