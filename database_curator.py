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
        line=line.rstrip()
        if '>' in line :
            seqs+=[seq]
            seq =''
        else:
            seq+=line
    return seqs
def reverse_complement(dna):
    dna=dna.replace('A','E').replace('T','A').replace('E','T').replace('G','S').replace('C','G').replace('S','C')
    return dna[::-1]
def remove_redudncy_N(totaldna,filename):          #remove redudncy returning the largest possible sequence
    f=open(filename,'w')
    while '' in totaldna:
        totaldna.remove('')
    comparingdna=''
    resultdna=[]
    totaldna.sort(key=len)
    for i in range(len (totaldna)):
        compileddna=''
        rcomparingdna=''
        comparingdna=totaldna[i]
        rcomparingdna=reverse_complement(comparingdna)
        compileddna=''.join(totaldna[(i+1):])
        if comparingdna in compileddna or rcomparingdna in compileddna:
            pass
        else:
            resultdna+=[comparingdna]
    for i in range(len(resultdna)):
        f.write(('\n>seq %d\n') %(i))
        f.write(resultdna[i])
    f.write('\n>')
    print 'deleted sequences =',(len(totaldna)+1-len(resultdna))
    f.close()
def remove_redudncy_Y(totaldna,filename,prot_length):             #remove redudncy returning the optimum length sequence
    f=open(filename,'w')
    while '' in totaldna:
        totaldna.remove('')
    comparingdna=''
    resultdna=[]
    totaldna.sort(key=len)
    for i in range(len (totaldna)):
        compileddna=''
        rcomparingdna=''
        comparingdna=totaldna[i]
        if comparingdna=='':
            pass
        elif (prot_length*2.7+2.7) <= len(comparingdna):
            rcomparingdna=reverse_complement(comparingdna)
            for j in range ((i+1),len (totaldna)):
                if comparingdna in totaldna[j] or rcomparingdna in totaldna[j]:
                    totaldna[j]=''
            resultdna+=[comparingdna]
        else:
            rcomparingdna=reverse_complement(comparingdna)
            compileddna=''.join(totaldna[(i+1):])
            if comparingdna in compileddna or rcomparingdna in compileddna:
                pass
            else:
                resultdna+=[comparingdna]
    for i in range(len(resultdna)):
        f.write(('\n>seq %d\n') %(i))
        f.write(resultdna[i])
    f.write('\n>')
    print 'deleted sequences =',(len(totaldna)+1-len(resultdna))
    f.close()
def remove_redudncy(total,filename):          #remove redudncy returning the largest possible sequence
    f=open(filename,'w')
    while '' in total:
        total.remove('')
    comparing=''
    result=[]
    total.sort(key=len)
    for i in range(len (total)):
        compiled=''
        comparing=total[i]
        compiled=''.join(total[(i+1):])
        if comparing in compiled:
            pass
        else:
            result+=[comparing]
    for i in range(len(result)):
        f.write(('\n>seq %d\n') %(i))
        f.write(result[i])
    f.write('\n>')
    print 'deleted sequences =',(len(total)+1-len(result))
    f.close()
database = str(raw_input("Protein or Nucleotide Database? P/N"))                  #begining of the code !!!
path = str(raw_input("Enter original FASTA file's path"))
file_name = str(raw_input("Enter original FASTA file's name"))
gene = str(raw_input("Enter gene's name"))
intermediate_file='%s_curated_seq_only.fasta' %gene
final_file='%s_final.fasta' %gene
import os
os.chdir(path)
import time                                                   #in case if we want the time
start_time = time.clock()
totalseqs=reading_FASTA_sequences(file_name)
medium_time=time.clock()                                        #in case if we want the time
if database.upper() == 'N':
    switch= str(raw_input("Do you want to have the optimum length approach? Y/N"))
    if switch.upper() == 'Y':
        prot_length = int(raw_input("Enter protein's length"))
        medium2_time=time.clock()                                        #in case if we want the time
        remove_redudncy_Y(totalseqs,intermediate_file,prot_length)
    elif switch.upper() == 'N':
        medium3_time=time.clock()                                        #in case if we want the time
        remove_redudncy_N(totalseqs,intermediate_file)
    else:
        print "Please specify your approach"
elif database.upper() == 'P':
    remove_redudncy(totalseqs,intermediate_file)
else:
    print "Please specify your type of database"
totalnames=reading_FASTA_names(file_name)
curatednames=[]
curatedseqs=reading_FASTA_sequences(intermediate_file)
del(curatedseqs[0])
f=open(final_file,'w')
for i in range(len(curatedseqs)):
    curatednames+=[totalnames[totalseqs.index(curatedseqs[i])]]
    f.write('\n%s\n%s' % (curatednames[i],curatedseqs[i]))
f.close()
if database=='P':
    time_of_calc = time.clock() - start_time
else:
    if switch.upper=='Y':
        time_of_calc= (time.clock()-medium2_time) + (medium_time-start_time)
    else:
        time_of_calc= (time.clock()-medium3_time) + (medium_time-start_time)
print time_of_calc, "seconds"                                      #in case if we want the time"""

