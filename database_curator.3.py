import sys
import os
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
    print ('deleted sequences =',(len(totalseqs)-len(result)))

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
    print ('deleted sequences =',(len(totalseqs)-len(result)))

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
    print ('deleted sequences =',(len(totalseqs)-len(result)))

                  #begining of the code !!!

database = str(input("Protein or Nucleotide Database? P/N : "))
if database.upper() != 'P' and database.upper() != 'N':
    sys.exit("\n\nPlease specify your type of database !\n")

path = str(input("Enter original FASTA file's path : "))
if path == '' or path == ' ':
    sys.exit("\n\nThere is no directory in this path !\n")

gene = str(input("Enter gene's name : "))
intermediate_file='%s_curated_seq_only.fasta' %gene
final_file='%s_final.fasta' %gene
deleted_file='%s_deleted.fasta' %gene
os.chdir(path)

multiples=str(input("Your data is in one file Y/N : "))
if multiples.upper() != 'Y' and multiples.upper() != 'N':
    sys.exit("\n\nPlease specify your data either in one file or not !\n")

if multiples.upper()=='N':
    ext=str(input("Enter your files' extension (fasta, fas, txt, ...) : "))
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
    file_name = str(input("Enter original FASTA file's name : "))
    if os.path.isfile(file_name) == False:
        sys.exit("\n\nYour file's name and / or your path is incorrect !\n")

start_time = time.clock()

totalnames,totalseqs=parsing(file_name)

medium_time=time.clock()                                        #in case if we want the time

if database.upper() == 'N':
    switch= str(input("Do you want to have the optimum length approach? Y/N : "))
    if switch.upper() == 'Y':
        prot_length = int(input("Enter protein's length : "))
        medium2_time=time.clock()                                        #in case if we want the time
        remove_redundancy_Y(totalseqs,totalnames,intermediate_file,prot_length)
    elif switch.upper() == 'N':
        medium3_time=time.clock()                                        #in case if we want the time
        remove_redundancy_N(totalseqs,totalnames,intermediate_file)
    else:
        sys.exit("\n\nPlease specify your approach !\n")
elif database.upper() == 'P':
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

if database.upper()=='P':
    time_of_calc = time.clock() - start_time
else:
    if switch.upper()=='Y':
        time_of_calc= (time.clock()-medium2_time) + (medium_time-start_time)
    else:
        time_of_calc= (time.clock()-medium3_time) + (medium_time-start_time)

print time_of_calc, "seconds"                                      #in case if we want the time"""
