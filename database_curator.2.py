import sys
import os
import argparse
import time
from Bio import SeqIO
from Bio.Seq import reverse_complement

def parsing (filename):                             #reading FASTA names and sequences
    records=SeqIO.parse(filename,'fasta')           #name,id,desciption,seq
    names,sequences=[],[]
    for record in records:
        names+=['>%s'%str(record.description)]
        sequences+=['%s'%str(record.seq)]
    return names, sequences

def remove_redundancy_N(totalseqs,totalnames,filename):          #remove redundancy returning the largest possible sequence
    editing=totalseqs[:]
    while '' in editing:
        editing.remove('')
    result=[]
    editing=list(set(editing))
    editing.sort(key=len)
    for i in range(len (editing)-1):
        comparing=editing[i]
        rcomparing=reverse_complement(comparing)
        compiled=''.join(editing[(i+1):])
        if comparing in compiled or rcomparing in compiled:
            pass
        else:
            result+=[comparing]
    result+=[editing[len (editing)-1]]
    with open(filename,'w') as f:
        f.write(('>seq %d\n%s') %((1),result[0]))
        for i in range(1,len(result)):
            f.write(('\n>seq %d\n%s') %((i+1),result[i]))
    print 'deleted sequences =',(len(totalseqs)-len(result))

def remove_redundancy_Y(totalseqs,totalnames,filename,prot_length):             #remove redundancy returning the optimum length sequence
    editing=totalseqs[:]
    while '' in editing:
        editing.remove('')
    result=[]
    editing=list(set(editing))
    editing.sort(key=len)
    for i in range(len (editing)):
        comparing=editing[i]
        if comparing=='':
            pass
        elif (prot_length*2.7+2.7) <= len(comparing):
            rcomparing=reverse_complement(comparing)
            for j in range ((i+1),len (editing)):
                if comparing in editing[j] or rcomparing in editing[j]:
                    editing[j]=''
            result+=[comparing]
        else:
            rcomparing=reverse_complement(comparing)
            compiled=''.join(editing[(i+1):])
            if comparing in compiled or rcomparing in compiled:
                pass
            else:
                result+=[comparing]
    with open(filename,'w') as f:
        f.write(('>seq %d\n%s') %((1),result[0]))
        for i in range(1,len(result)):
            f.write(('\n>seq %d\n%s') %((i+1),result[i]))
    print 'deleted sequences =',(len(totalseqs)-len(result))

def remove_redundancy(totalseqs,totalnames,filename):          #remove redundancy returning the largest possible sequence
    editing=totalseqs[:]
    while '' in editing:
        editing.remove('')
    result=[]
    editing=list(set(editing))
    editing.sort(key=len)
    for i in range(len (editing)-1):
        comparing=editing[i]
        compiled=''.join(editing[(i+1):])
        if comparing in compiled:
            pass
        else:
            result+=[comparing]
    result+=[editing[len (editing)-1]]
    with open(filename,'w') as f:
        f.write(('>seq %d\n%s') %((1),result[0]))
        for i in range(1,len(result)):
            f.write(('\n>seq %d\n%s') %((i+1),result[i]))
    print 'deleted sequences =',(len(totalseqs)-len(result))

                  #begining of the code !!!

parser = argparse.ArgumentParser(prog='Database curator program',usage='\n%(prog)s :curates nucleotide and/or protein databases from redundant and partial redundant sequences.\n\n Eslam S.Ibrahim\n\n eslam.ebrahim@pharma.cu.edu.eg')
parser.add_argument('-p',dest='database',action='store_const',const='p',help='protein sequences')
parser.add_argument('-n',dest='database',action='store_const',const='n',help='nucleotide sequences')
parser.add_argument('-desired',type=str,dest='gene',help='desired name for your files')
parser.add_argument('-multi',dest='multiples',default=False,action='store_true',help='if there are multiple files to process')
parser.add_argument('-in',dest='file_name',type=argparse.FileType('r'),required=True)
parser.add_argument('-optimum',dest='optimum_length_approach',default=False,action='store_true',help='if optimum length approach are wanted')
args=parser.parse_args()

database = args.database
gene = args.gene
multiples=args.multiples

if database != 'p' and database != 'n':
    sys.exit("\n\nPlease specify your type of database !\n")

intermediate_file='%s_curated_seq_only.fasta' %gene
final_file='%s_final.fasta' %gene
deleted_file='%s_deleted.fasta' %gene

if multiples==True:
    path = str(raw_input("Enter original FASTA files' path : "))
    ext=str(raw_input("Enter your files' extension (fasta, fas, txt, ...) : "))
    if '.' in ext:
        sys.exit("\n\nPlease write the extension only without '.' !\n")
    import glob
    files = glob.glob('*.%s'%ext)
    if len(files) == 0:
        sys.exit("\n\nYour files' extension and / or your path is incorrect !\n")
    with open('%s.%s'%(gene,ext),'w') as result:
        for _ in files:
            for line in open(_,'r'):
                result.write(line)
    file_name='%s.fasta' %gene
else:
    file_name = args.file_name

start_time = time.clock()

totalnames,totalseqs=parsing(file_name)

medium_time=time.clock()                                        #in case if we want the time

if database == 'n':
    optimum = args.optimum_length_approach
    if optimum == True:
        prot_length = int(raw_input("Enter protein's length : "))
        medium2_time=time.clock()                                        #in case if we want the time
        remove_redundancy_Y(totalseqs,totalnames,intermediate_file,prot_length)
    elif optimum == False:
        medium3_time=time.clock()                                        #in case if we want the time
        remove_redundancy_N(totalseqs,totalnames,intermediate_file)
    else:
        sys.exit("\n\nPlease specify your approach !\n")
elif database == 'p':
    remove_redundancy(totalseqs,totalnames,intermediate_file)

curatedseqs=parsing(intermediate_file)[1]
curatednames=[totalnames[totalseqs.index(seq)] for seq in curatedseqs]
deleted= list(set([name for name in totalnames if name not in curatednames]))

with open(final_file,'w') as f:
    f.write('%s\n%s' % (curatednames[0],curatedseqs[0]))
    for i in range(1,len(curatedseqs)):
        f.write('\n%s\n%s' % (curatednames[i],curatedseqs[i]))

with open(deleted_file,'w') as f:
    f.write('%s\n'%deleted[0])
    for i in range(1,len(deleted)):
        f.write('\n%s\n'%deleted[i])

if database=='p':
    time_of_calc = time.clock() - start_time
else:
    if optimum==True:
        time_of_calc= (time.clock()-medium2_time) + (medium_time-start_time)
    else:
        time_of_calc= (time.clock()-medium3_time) + (medium_time-start_time)

print time_of_calc, "seconds"                                      #in case if we want the time"""
