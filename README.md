# database_curator
This program can curate nucleotide and/or protein databases from redundant and partial redundant sequences for a specific gene and /or any groups of genes.

# Input:
File containing all the different downloaded sequences in FASTA format.

# Processing:
1- It removes the redundant sequences.
2- It removes the partial sequences that are exact part from other sequences in your database.

# Output:
Two files (%gene_curated_seq_only.fasta and %gene_final.fasta) contain only the unique and more complete sequences in your input file depending on your approach.
%gene_curated_seq_only.fasta : contains sequences without their annotations.
%gene_final.fasta : contains sequences with the exact annotation in your database

# Options:
1- Working on either Protein (P) or Nucleotide (N) databases.
2- Two approaches (largest possible length and optimum length).
   2a- largest possible length approach: gives the longest sequence even if it exceeds the length of your gene.
   2b- optimum length approach: gives only your gene provided you feed the approximate length of your protein.

# How to use:
1-	you need to install python 2.7 (https://www.python.org/downloads/) on your machine.
2-	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
3-	Open the file “database_curator.py” with (python.exe).
4-	State your variables and press Enter.


# Any Errors please let me know via my e-mail: eslam.ebrahim@pharma.cu.edu.eg
