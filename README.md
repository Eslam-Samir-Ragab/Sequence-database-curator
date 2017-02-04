# Sequence database curator
This program can curate nucleotide and/or protein databases from redundant and partial redundant sequences for a specific gene and /or any groups of genes.
![Summary of the program.](https://ppjpkw-ch3302.files.1drv.com/y3pdwEcVS3Zzm79n26oygTt4nJ1R-52Ziqvxpr29UPFhd72AERdJvinl95KnTvT5WZ8m9iScsinX0cLTyyMdly24jYX1sebuuEzzFuGpTlkxGwbtpwlU_LHZDJrVT3gHSV6ub79tHDQXg-d_tNT1GVWwqVbZ67YP-QDCosfs2zBuUo/Database%20curator.png?psid=1)

## Input:
- File containing all the different downloaded sequences in FASTA format.
- Or a directory containing all your desired files in the same extension.

## Processing:
1. It removes the redundant sequences.
2. It removes the partial sequences that are exact part from other sequences in your database.
## Output:
Three files (**%gene_curated_seq_only.fasta**, **%gene_final.fasta** and **%gene_deleted.fasta**) contain only the unique and more complete sequences in your input file depending on your approach.
- *%gene_curated_seq_only.fasta* : contains sequences without their annotations.
- *%gene_final.fasta* : contains sequences with the exact annotation in your database.
- *%gene_deleted.fasta* : contains the names of the deleted seqeunces.

## Options:
1. Working on either **Protein (P) or Nucleotide (N)** databases.
2. Two approaches (**largest possible length** and **optimum length**).
  * largest possible length approach: gives the longest sequence even if it exceeds the length of your gene.
  * optimum length approach: gives only your gene provided you feed the approximate length of your protein.
![Illustration of both approaches.](https://ppjqaa-ch3302.files.1drv.com/y3p6MyxtnFjVwWErixgUKwFIo5p2TQTrMCdzkWUTBK8yDPWhyqeTJHC8bZwrO1dx1PE9Whj6pKaPSpWg3eiUSNhM59AZBre77KnE7QS95ME1MP7GSne3DjOlJo_0e2JgR_JPLNgR69UHSoZxNPjs0ZY7qEO6utxPfU93PFp7uxMubI/Capture%20%281%29.PNG?psid=1)

## How to use (Python 2.7 or later):
1.	you need to install [python 2.7](https://www.python.org/downloads/) on your machine.
2.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
3.	Open the file “**database_curator.2.py**” with (python.exe).
  * [Windows](http://stackoverflow.com/a/1527012/7414020)
  * U/Linux : use the file “**database_curator.linux.2.py**” and the command `chmod u+x database_curator.linux.2.py`
  * Mac : use the file “**database_curator.2.py**” and the command `python database_curator.2.py`
4.	State your variables and press Enter.

## How to use (Python 3):
1.	you need to install [python 3](https://www.python.org/downloads/) on your machine.
2.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
3.	Open the file “**database_curator.3.py**” with (python.exe).
  * **[Windows](http://stackoverflow.com/a/1527012/7414020)**
  * **U/Linux** : use the file “**database_curator.linux.3.py**” and the command `chmod u+x database_curator.linux.3.py`
  * **Mac** : use the file “**database_curator.2.py**” and the command `python database_curator.3.py`
4.	State your variables and press Enter.


### Any errors please let me know via an e-mail with the subject "database_curator" to eslam.ebrahim@pharma.cu.edu.eg

# Fastq Dereplicator
This program can remove Fastq sequences from one or multiple fastq files and return the dereplicated version in fasta format for other processing.
## Input:
- File containing all the different sequences in fastq format.
- Or a directory containing all your desired files in the same extension.

## Processing:
- It removes the redundant sequences.

## Output:
Three files (**?_curated_seq_only.fasta**, **?_final.fasta** and **?.fasta**) contain only the unique sequences in your input file depending on your approach.
- *?_curated_seq_only.fasta* : contains sequences without their annotations.
- *?_final.fasta* : contains sequences with the exact annotation in your input file.
- *?.fasta* : the same your fastq file but in a FASTA format.

## Options:
1. Two approaches (**largest possible length** and **smallest possible length**).
  * largest possible length approach: gives the longest sequence.
  * smallest possible length: gives the shortest sequence in your dataset.

## How to use (Python 2.7 or later):
1.	you need to install [python 2.7](https://www.python.org/downloads/) on your machine.
2. you need to install [Numpy](https://pypi.python.org/pypi/numpy) and [Biopython](http://biopython.org/wiki/Download)
3.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
4.	Open the file “**fastq_dereplicator.py**” with (python.exe).
  * [Windows](http://stackoverflow.com/a/1527012/7414020)
  * U/Linux : use the file “**fastq_dereplicator.linux.py**” and the command `chmod u+x fastq_dereplicator.py`
  * Mac : use the file “**fastq_dereplicator.py**” and the command `python fastq_dereplicator.py`
4.	State your variables and press Enter.
