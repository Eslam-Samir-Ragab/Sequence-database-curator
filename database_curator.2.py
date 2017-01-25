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
def remove_redudncy_N(totalseqs,totalnames,filename,filename2):          #remove redudncy returning the largest possible sequence
    f=open(filename,'w')
    editing=totalseqs[:]
    while '' in editing:
        editing.remove('')
    result=[]
    deleted=[]
    editing.sort(key=len)
    for i in range(len (editing)-1):
        comparing=editing[i]
        rcomparing=reverse_complement(comparing)
        compiled=''.join(editing[(i+1):])
        if comparing in compiled or rcomparing in compiled:
            deleted+=[totalnames[totalseqs.index(comparing)]]
        else:
            result+=[comparing]
    result+=[editing[len (editing)-1]]
    f.write(('>seq %d\n%s') %((1),result[0]))
    for i in range(1,len(result)):
        f.write(('\n>seq %d\n%s') %((i+1),result[i]))
    m=open(filename2,'w')
    m.write('%s\n'%deleted[0])
    for i in range(1,len(deleted)):
        m.write('\n%s\n'%deleted[i])
    m.close()
    print 'deleted sequences =',(len(totalseqs)-len(result))
    f.close()
def remove_redudncy_Y(totalseqs,totalnames,filename,filename2,prot_length):             #remove redudncy returning the optimum length sequence
    f=open(filename,'w')
    editing=totalseqs[:]
    while '' in editing:
        editing.remove('')
    result=[]
    deleted=[]
    editing.sort(key=len)
    for i in range(len (editing)):
        comparing=editing[i]
        if comparing=='':
            pass
        elif (prot_length*2.7+2.7) <= len(comparing):
            rcomparing=reverse_complement(comparing)
            for j in range ((i+1),len (editing)):
                if comparing in editing[j] or rcomparing in editing[j]:
                    deleted+=[totalnames[totalseqs.index(editing[j])]]
                    editing[j]=''
            result+=[comparing]
        else:
            rcomparing=reverse_complement(comparing)
            compiled=''.join(editing[(i+1):])
            if comparing in compiled or rcomparing in compiled:
                deleted+=[totalnames[totalseqs.index(comparing)]]
            else:
                result+=[comparing]
    f.write(('>seq %d\n%s') %((1),result[0]))
    for i in range(1,len(result)):
        f.write(('\n>seq %d\n%s') %((i+1),result[i]))
    m=open(filename2,'w')
    m.write('%s\n'%deleted[0])
    for i in range(1,len(deleted)):
        m.write('\n%s\n'%deleted[i])
    m.close()
    print 'deleted sequences =',(len(totalseqs)-len(result))
    f.close()
def remove_redudncy(totalseqs,totalnames,filename,filename2):          #remove redudncy returning the largest possible sequence
    f=open(filename,'w')
    editing=totalseqs[:]
    while '' in editing:
        editing.remove('')
    result=[]
    deleted=[]
    editing.sort(key=len)
    for i in range(len (editing)-1):
        comparing=editing[i]
        compiled=''.join(editing[(i+1):])
        if comparing in compiled:
            deleted+=[totalnames[totalseqs.index(comparing)]]
        else:
            result+=[comparing]
    result+=[editing[len (editing)-1]]
    f.write(('>seq %d\n%s') %((1),result[0]))
    for i in range(1,len(result)):
        f.write(('\n>seq %d\n%s') %((i+1),result[i]))
    m=open(filename2,'w')
    m.write('%s\n'%deleted[0])
    for i in range(1,len(deleted)):
        m.write('\n%s\n'%deleted[i])
    m.close()
    print 'deleted sequences =',(len(totalseqs)-len(result))
    f.close()
database = str(raw_input("Protein or Nucleotide Database? P/N : "))                  #begining of the code !!!
path = str(raw_input("Enter original FASTA file's path : "))
gene = str(raw_input("Enter gene's name : "))
intermediate_file='%s_curated_seq_only.fasta' %gene
final_file='%s_final.fasta' %gene
deleted_file='%s_deleted.fasta' %gene
import os
os.chdir(path)
multiples=str(raw_input("Your data is in one file Y/N : "))
if multiples.upper()=='N':
    ext=str(raw_input("Enter your files' extention (fasta, fas, txt, ...) : "))
    import glob
    files = glob.glob('*.%s'%ext)
    with open('%s.%s'%(gene,ext),'w') as result:
        for _ in files:
            for line in open(_,'r'):
                result.write(line)
    file_name='%s.fasta' %gene
else:
    file_name = str(raw_input("Enter original FASTA file's name : "))
import time                                                   #in case if we want the time
start_time = time.clock()
totalseqs=reading_FASTA_sequences(file_name)
medium_time=time.clock()                                        #in case if we want the time
totalnames=reading_FASTA_names(file_name)
if database.upper() == 'N':
    switch= str(raw_input("Do you want to have the optimum length approach? Y/N : "))
    if switch.upper() == 'Y':
        prot_length = int(raw_input("Enter protein's length : "))
        medium2_time=time.clock()                                        #in case if we want the time
        remove_redudncy_Y(totalseqs,totalnames,intermediate_file,deleted_file,prot_length)
    elif switch.upper() == 'N':
        medium3_time=time.clock()                                        #in case if we want the time
        remove_redudncy_N(totalseqs,totalnames,intermediate_file,deleted_file)
    else:
        print "Please specify your approach"
elif database.upper() == 'P':
    remove_redudncy(totalseqs,totalnames,intermediate_file,deleted_file)
else:
    print "Please specify your type of database"
curatednames=[]
curatedseqs=reading_FASTA_sequences(intermediate_file)
f=open(final_file,'w')
for i in range(len(curatedseqs)):
    curatednames+=[totalnames[totalseqs.index(curatedseqs[i])]]
f.write('%s\n%s' % (curatednames[0],curatedseqs[0]))
for i in range(1,len(curatedseqs)):
    f.write('\n%s\n%s' % (curatednames[i],curatedseqs[i]))
f.close()
if database.upper()=='P':
    time_of_calc = time.clock() - start_time
else:
    if switch.upper()=='Y':
        time_of_calc= (time.clock()-medium2_time) + (medium_time-start_time)
    else:
        time_of_calc= (time.clock()-medium3_time) + (medium_time-start_time)
print time_of_calc, "seconds"                                      #in case if we want the time"""
