# Sequence database curator
This program can curate nucleotide and/or protein databases from redundant and partial redundant sequences for a specific gene and /or any groups of genes.

## Input:
File containing all the different downloaded sequences in FASTA format.

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
![Illustration of both approaches.](https://drive.google.com/file/d/0B752ugm5kNQcTGJGd203blk5UTQ/view?usp=sharing)

## How to use (Python 2.7 or later):
1.	you need to install [python 2.7](https://www.python.org/downloads/) on your machine.
2.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
3.	Open the file “**database_curator.2.py**” with (python.exe).
  * [Windows](http://stackoverflow.com/a/1527012/7414020)
  * U/Linux : use the file “**(Unix-Linux)database_curator.2.py**” and the command `chmod u+x (Unix-Linux)database_curator.2.py`
  * Mac : use the file “**database_curator.2.py**” and the command `python database_curator.2.py`
4.	State your variables and press Enter.

## How to use (Python 3):
1.	you need to install [python 3](https://www.python.org/downloads/) on your machine.
2.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
3.	Open the file “**database_curator.3.py**” with (python.exe).
  * **[Windows](http://stackoverflow.com/a/1527012/7414020)**
  * **U/Linux** : use the file “**(Unix-Linux)database_curator.3.py**” and the command `chmod u+x (Unix-Linux)database_curator.3.py`
  * **Mac** : use the file “**database_curator.2.py**” and the command `python database_curator.3.py`
4.	State your variables and press Enter.


### Any errors please let me know via an e-mail with the subject "database_curator" to eslam.ebrahim@pharma.cu.edu.eg
