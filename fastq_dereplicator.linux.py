#!/usr/bin/env python
def reading_FASTA_names(filename):                 #reading FASTA names
    f=open(filename)
    all_lines=f.readlines()
    f.close()
    names=[all_lines[0].rstrip()]
    for i in range (1,len (all_lines)):
        line=all_lines[i]
        line=line.rstrip()
        if '>' in line :
            names.append(line)
    return names
def reading_FASTA_sequences(filename):              #reading FASTA sequences
    f=open(filename)
    all_lines=f.readlines()
    f.close()
    seqs=[]
    seq=''
    for i in range (1,len (all_lines)):
        line=all_lines[i]
        line=line.replace(' ','')
        line=line.rstrip()
        if '>' in line :
            seqs+=[seq]
            seq =''
        else:
            seq+=line
    seqs+=[seq]
    return seqs
def reverse_complement(dna):
    dna=dna.replace('A','E').replace('T','A').replace('E','T').replace('G','S').replace('C','G').replace('S','C')
    return dna[::-1]
def remove_redudncy_N(totalseqs,totalnames,filename):          #remove redudncy returning the largest possible sequence
    f=open(filename,'w')
    editing=totalseqs[:]
    while '' in editing:
        editing.remove('')
    result=[]
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
    f.write(('>seq %d\n%s') %((1),result[0]))
    for i in range(1,len(result)):
        f.write(('\n>seq %d\n%s') %((i+1),result[i]))
    print 'deleted sequences =',(len(totalseqs)-len(result))
    f.close()
def remove_redudncy_Y(totalseqs,totalnames,filename):             #remove redudncy returning the smallest length sequence
    f=open(filename,'w')
    editing=totalseqs[:]
    while '' in editing:
        editing.remove('')
    result=[]
    editing.sort(key=len)
    for i in range(len (editing)):
        comparing=editing[i]
        if comparing=='':
            pass
        else:
            rcomparing=reverse_complement(comparing)
            for j in range ((i+1),len (editing)):
                if comparing in editing[j] or rcomparing in editing[j]:
                    editing[j]=''
            result+=[comparing]
    f.write(('>seq %d\n%s') %((1),result[0]))
    for i in range(1,len(result)):
        f.write(('\n>seq %d\n%s') %((i+1),result[i]))
    print 'deleted sequences =',(len(totalseqs)-len(result))
    f.close()

                 #begining of the code !!!

path = str(raw_input("Enter original Fastq file's path : "))
if path == '' or path == ' ':
    import sys
    sys.exit("\n\nThere is no directory in this path !\n")
desired = str(raw_input("Enter your desired files' name : "))
intermediate_file='%s_curated_seq_only.fasta' %desired
final_file='%s_final.fasta' %desired
deleted_file='%s_deleted.fasta' %desired
import os
os.chdir(path)
multiples=str(raw_input("Your data is in one file Y/N : "))
if multiples.upper() != 'Y' and multiples.upper() != 'N':
    import sys
    sys.exit("\n\nPlease specify your data either in one file or not !\n")
if multiples.upper()=='N':
    import glob
    files = glob.glob('*.fastq')
    if len(files) == 0:
        import sys
        sys.exit("\n\nYour files' extension and / or your path is incorrect !\n")
    with open('%s.fastq'%(desired),'w') as result:
        for _ in files:
            for line in open(_,'r'):
                result.write(line)
    fastq_file_name='%s.fastq' %desired
else:
    fastq_file_name = str(raw_input("Enter original Fastq file's name : "))
    if os.path.isfile(fastq_file_name) == False:
        import sys
        sys.exit("\n\nYour file's name and / or your path is incorrect !\n")
from Bio import SeqIO
import time                                                   #in case if we want the time
start_time = time.clock()
file_name="%s.fasta"%desired
SeqIO.convert(fastq_file_name, "fastq", file_name, "fasta")
totalseqs=reading_FASTA_sequences(file_name)
medium_time=time.clock()                                        #in case if we want the time
totalnames=reading_FASTA_names(file_name)
switch= str(raw_input("Do you want to have the smallest length approach? Y/N : "))
if switch.upper() == 'Y':
    medium2_time=time.clock()                                        #in case if we want the time
    remove_redudncy_Y(totalseqs,totalnames,intermediate_file)
elif switch.upper() == 'N':
    medium3_time=time.clock()                                        #in case if we want the time
    remove_redudncy_N(totalseqs,totalnames,intermediate_file)
else:
    import sys
    sys.exit("\n\nPlease specify your approach !\n")
curatednames=[]
curatedseqs=reading_FASTA_sequences(intermediate_file)
f=open(final_file,'w')
for i in range(len(curatedseqs)):
    curatednames+=[totalnames[totalseqs.index(curatedseqs[i])]]
f.write('%s\n%s' % (curatednames[0],curatedseqs[0]))
for i in range(1,len(curatedseqs)):
    f.write('\n%s\n%s' % (curatednames[i],curatedseqs[i]))
f.close()
if switch.upper()=='Y':
    time_of_calc= (time.clock()-medium2_time) + (medium_time-start_time)
else:
    time_of_calc= (time.clock()-medium3_time) + (medium_time-start_time)
print time_of_calc, "seconds"                                      #in case if we want the time"""
