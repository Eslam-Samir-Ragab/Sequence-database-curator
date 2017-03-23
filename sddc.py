#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import argparse
import time
from Bio import SeqIO
from Bio.Seq import reverse_complement

def parsing_file (filename):                             #reading FASTA names and sequences
    records=SeqIO.parse(filename,'fasta')
    names,sequences=[],[]
    for record in records:
        names+=['>%s'%str(record.description)]
        sequences+=['%s'%str(record.seq)]
    return names, sequences

def file_writer(output_file,starting_names,starting_sequences,curated_sequences):           #write function
    print ('\n-------\nfiltered sequences = %d from %d starting sequences\nresulting sequence = %d sequences' %(len(starting_names)-len(curated_sequences),len(starting_names),len(curated_sequences)))
    with open(output_file,'w') as f:
        name=starting_names[starting_sequences.index(curated_sequences[0])]
        curated_names=[name]
        f.write('%s\n%s' % (name,curated_sequences[0]))
        for i in range(1,len(curated_sequences)):
            name=starting_names[starting_sequences.index(curated_sequences[i])]
            curated_names.append(name)
            f.write('\n%s\n%s' % (name,curated_sequences[i]))
    deleted= list(set([name for name in starting_names if name not in curated_names]))
    if len(deleted) > 0:
        with open('names_of_deleted.txt','w') as f:
            f.write('%s\n'%deleted[0])
            for i in range(1,len(deleted)):
                f.write('\n%s\n'%deleted[i])

def cleaner(_list):                                      #to clean after processing
    while '' in _list:
        _list.remove('')
    return _list

def kmer_gen(sequence,k,start):                #to generate k-mers
    kmer=sequence[start:start+k+1]
    return kmer

def remover_by_seq(start_file,remove_file,database):       #remover (by sequences)
    starting_names,starting_sequences=parsing_file(start_file)
    removing_sequences=set(parsing_file(remove_file)[1])
    if len(starting_sequences) < 1 or len(removing_sequences) < 1:
        sys.exit("\n\nError in your input file or filteration file or incorrect filteration mode !\n")

    editing = set (starting_sequences[:])
    common = editing & removing_sequences
    removing_sequences = list(removing_sequences-common)
    editing = list(editing-common)

    for i in range(len(editing)):
        comparing=editing[i]
        for seq in removing_sequences:
            kmer = kmer_gen(comparing,length,0)
            if kmer in seq:
                if comparing in seq:
                    editing[i]=''
            elif database=='n' and reverse_complement(kmer) in seq:
                if reverse_complement(comparing) in seq:
                    editing[i]=''
    editing=cleaner(editing)
    file_writer(output_file,starting_names,starting_sequences,editing)

def remover_by_name (start_file,remove_file,filteration):     #remover (by names)
    with open(remove_file,'r') as f:
        lines=f.readlines()
        filter_names=[line.rstrip() for line in lines if '>' in line]

    starting_names,starting_sequences=parsing_file(start_file)
    data=[(starting_names[i]+'&&'+starting_sequences[i]) for i in range(len(starting_names))]
    data=list(set(data))
    if filteration == 'exclusive':
        curated_data=[line for line in data for filter_name in filter_names if filter_name not in line]
        curated_data=list(set(curated_data))
    elif filteration =='inclusive':
        curated_data=[line for line in data for filter_name in filter_names if filter_name in line]
    print ('\n-------\nfiltered sequences = %d from %d starting sequences\nresulting sequence = %d sequences' %(len(starting_names)-len(curated_data),len(starting_names),len(curated_data)))
    with open(output_file,'w') as f:
        lines_to_write = [line.split('&&') for line in curated_data]
        for i in range(len(lines_to_write)):
            f.write('\n%s\n%s' % (lines_to_write[i][0],lines_to_write[i][1]))

def derep_longest(start_file,database):                         #dereplication (longest approach)
    starting_names,starting_sequences=parsing_file(start_file)
    if len(starting_sequences) < 1 :
        sys.exit("\n\nError in your input file !\n")

    editing=list(set(starting_sequences))
    editing.sort(key=len)

    for i in range((len(editing)-1)):
        if len(editing[i]) >= minimum:
            comparing=editing[i]
            kmer = kmer_gen(comparing,length,0)
            for seq in editing[(i+1):]:
                if kmer in seq:
                    if comparing in seq:
                        editing[i]=''
                elif database=='n' and reverse_complement(kmer) in seq:
                    if reverse_complement(comparing) in seq:
                        editing[i]=''

    editing=cleaner(editing)
    file_writer(output_file,starting_names,starting_sequences,editing)

def derep_optimum(start_file,prot_length):                     #dereplication (optimum approach)
    starting_names,starting_sequences=parsing_file(start_file)
    if len(starting_sequences) < 1:
        sys.exit("\n\nError in your input file !\n")

    editing=list(set(starting_sequences))
    editing=cleaner(editing)
    editing.sort(key=len)

    for i in range(len(editing)):
        comparing= editing[i]
        for j in range((i+1),len (editing)):
            seq=editing[j]
            if len(editing[i]) >= minimum:
                kmer=kmer_gen(comparing,length,0)
                if kmer in seq:
                    if comparing in seq and (prot_length*2.7+2.7) >= len(seq):
                        editing[i]=''
                        break
                    elif comparing in seq and (prot_length*2.7+2.7) <= len(seq):
                        editing[j]=''
                if database =='n' and reverse_complement(kmer) in seq:
                    if reverse_complement(comparing) in seq and (prot_length*2.7+2.7) >= len(seq):
                        editing[i]=''
                        break
                    elif reverse_complement(comparing) in seq and (prot_length*2.7+2.7) <= len(seq):
                        editing[j]=''
            else:
                editing[i]=''
                break
    editing=cleaner(editing)
    file_writer(output_file,starting_names,starting_sequences,editing)


                  #begining of the code !!!

parser = argparse.ArgumentParser(prog='Sequence Database Dereplicator and Curator (SDDC) program',usage='\n%(prog)s : dereplicates and/or filter nucleotide and/or protein database from a list of names or sequences (by exact match).\n\n Eslam S.Ibrahim\n\n eslam.ebrahim@pharma.cu.edu.eg')
parser.add_argument('-mode',dest='mode',required=True,choices=['derep','filter'],help='dereplicate your file/files or filter your file from specific sequences or names')
parser.add_argument('-approach',dest='filter',choices=['inclusive','exclusive'],default='exclusive',help='if you want to filter your file/files by names and/or sequences either inclusively or exclusively')
parser.add_argument('-out',dest='output_file',required=True,help='Your output file')
parser.add_argument('-in',dest='input_file',type=argparse.FileType('r'),required=True,help='Input file containing your original data to be dereplicated and/or filtered')
parser.add_argument('-flt_file',dest='flt_file',help='Input file containing your listed names or sequences to be filtered from your original file')
parser.add_argument('-flt_by',dest='filter_approach',choices=['seq','name'],default='seq',help='The approach by which you wan your input file to be filtered')
parser.add_argument('-p',dest='database',action='store_const',const='p',help='protein sequences')
parser.add_argument('-n',dest='database',action='store_const',const='n',help='nucleotide sequences')
parser.add_argument('-len',dest='length',default=20,type=int,help='length of k-mer (for seq filter mode only)')
parser.add_argument('-min_length',dest='minimum',default=1,type=int,help='minimum sequence length in your data (for dereplication only)')
parser.add_argument('-multi',dest='multiples',default=False,action='store_true',help='if there are multiple files to process')
parser.add_argument('-fastq',dest='fastq',default=False,action='store_true',help='if your data is in fastq format')
parser.add_argument('-optimum',dest='optimum_length_approach',default=False,action='store_true',help='if optimum length approach is wanted')
parser.add_argument('-prot_length',dest='prot_length',type=int,help='protein length (for optimum approach in dereplication mode only)')
args=parser.parse_args()

                  #Parsing the argument !!!

mode=args.mode
filteration=args.filter
output_file=args.output_file
input_file=args.input_file
remove_file=args.flt_file
filter_approach=args.filter_approach
database=args.database
length=args.length
minimum=args.minimum
multiples=args.multiples
fastq=args.fastq
optimum = args.optimum_length_approach
prot_length = args.prot_length

                  #Checkers !!!

start_time = time.clock()

if database != 'p' and database != 'n':
    sys.exit("\n\nPlease specify your type of database !\n")

if (filteration == 'inclusive' or filteration == 'exclusive') and remove_file=='':
    sys.exit("\n\nError with filteration file or filteration criteria (inclusive or exclusive)!\n")

if filteration == 'inclusive' and filter_approach=='seq':
    sys.exit("\n\nfilteration criteria (inclusive or exclusive) works only if filter mode is by names\n")

if length < 20:
    sys.exit("\n\nminimum k-mer length = 20\n")

if optimum and prot_length < 1:
    sys.exit("\n\nPlease enter your portein length by the flag (-prot_length) !\n")


                  #Processing !!!

if mode=='filter' and filter_approach=='name':
    remover_by_name(input_file,remove_file,filteration)

elif mode=='filter' and filter_approach=='seq':
    remover_by_seq(input_file,remove_file,database)

elif mode=='derep' and multiples:
    path = str(input("Enter original FASTA files' path : "))
    ext=str(input("Enter your files' extension (fasta, fas, txt, ...) : "))
    if '.' in ext:
        sys.exit("\n\nPlease write the extension only without '.' !\n")
    import glob
    files = glob.glob('*.%s'%ext)
    if len(files) == 0:
        sys.exit("\n\nYour files' extension and / or your path is incorrect !\n")
    with open('compiled.%s'%ext,'w') as f:
        for _ in files:
            for line in open(_,'r'):
                f.write(line)
    input_file='compiled.%s'%ext
    if optimum:
        derep_optimum(input_file,prot_length)
    else:
        derep_longest(input_file,database)

elif mode=='derep' and fastq:
    SeqIO.convert(input_file, "fastq", 'converted.fasta', "fasta")
    input_file='converted.fasta'
    if optimum:
        derep_optimum(input_file,prot_length)
    else:
        derep_longest(input_file,database)

elif mode=='derep' and optimum:
    derep_optimum(input_file,prot_length)

elif mode=='derep' and optimum==False:
    derep_longest(input_file,database)

                  #Finishing !!!

time_of_calc = time.clock() - start_time
print (time_of_calc, "seconds")                                      #in case if we want the time
print ('-------\nThanks for using SDDC\n-------\nfor contact ==> eslam.ebrahim@pharma.cu.edu.eg\n-------')
