# database_curator
This program can curate nucleotide and/or protein databases from redundant and partial redundant seqeunces for a specific gene and /or any groups of genes

1- It accepts FASTA formats of both types of databases (in one file).
2- It removes the redundent seqeunces.
3- It removes the partial seqeunces that are exact part from other sequences in your database, so you have only the most complete seqeunce in your database.
4- It gives you two files (%gene_curated_seq_only.fasta and gene_final.fasta).
5- The first file is the sequences without it's name.
6- The second file is the seqeunces with the exact annotaion in your database.
7- It gives you the number of deleted seqeunces and the time of running the program.

Options:
1- Working on both Protein (P) or Nucleotide (N) databases.
2- Two approaches (largest possible length and optimum length).
   2a- largest possible length approach : gives the longest sequence even if it exceeds the length of your gene.
   2b- optimum length approach : gives only your gene provided you feed the approximate length of your protein.
