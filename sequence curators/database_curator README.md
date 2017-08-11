# Sequence database curator
This program can curate nucleotide and/or protein databases from redundant and partial redundant sequences for a specific gene and /or any groups of genes.

## Please, cite: DOI: 10.13140/RG.2.2.24200.11529


<a href="http://www.freeimagehosting.net/commercial-photography/"><img src="http://i.imgur.com/eFhxVTF.png" alt="commercial photography locations"></a>

## Input:
- File containing all the different downloaded sequences in FASTA or FASTQ format.
- Or a directory containing all your desired files in the same extension.

## Processing:
1. It removes the redundant sequences.
2. It removes the partial sequences that are exact part from other sequences in your database.
## Output:
Three files (**curated_seq_only.fasta**, **final.fasta** and **deleted.fasta**) contain only the unique and more complete sequences in your input file depending on your approach.
- *curated_seq_only.fasta* : contains sequences without their annotations.
- *final.fasta* : contains sequences with the exact annotation in your database.
- *deleted.fasta* : contains the names of the deleted seqeunces.
- *.fasta* : contains your data in fasta format if your data is in fastq.

## Options:
1. Working on either **protein (p) or nucleotide (n)** databases.
2. Two approaches (**largest possible length** and **optimum length**).
  * largest possible length approach: gives the longest sequence even if it exceeds the length of your gene.
  * optimum length approach: gives only your gene provided you feed the approximate length of your protein.
<a href="http://www.freeimagehosting.net/commercial-photography/"><img src="http://i.imgur.com/H0EOUf8.png" alt="commercial photography locations"></a>

## How to use:
1.	you need to install [python 2.7](https://www.python.org/downloads/) or [python 3](https://www.python.org/downloads/) on your machine.
2. you need to install [Numpy](https://pypi.python.org/pypi/numpy) and [Biopython](http://biopython.org/wiki/Download)
3. you need to install future module by [pip command](https://docs.python.org/3/installing/)
4.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
5.	Open the file “**database_curator.py**” with (python.exe).
  * [Windows](http://stackoverflow.com/a/1527012/7414020)
  * U/Linux : use the command `chmod u+x database_curator.py`
  * Mac : use the command `python database_curator.py`
6.	State your variables and press Enter.

List of options in the program:

| No |    Option   |                                            function                                           |
|:--:|:-----------:|:---------------------------------------------------------------------------------------------:|
|  1 |   * -in     | indicate input file, you can use the full path to your file (C:\Users\user\Desktop\dna.fasta) |
|  2 |   *  -n     | process nucleotide sequences                                                                  |
|  3 |   *  -p     | process protein sequences                                                                     |
|  4 |   -desired  | indicate a nomenclature for your output files                                                 |
|  5 |   -fastq    | if your data is in fastq format                                                               |
|  6 |    -multi   | if there are multiple files to process                                                        |
|  7 |   -optimum  | if optimum length approach is wanted'                                                         |
|  8 | -h or -help | show the help for the program                                                                 |

**Required** (*)  -n and -p are mutually exclusive flags


## Examples

if you want to process a nucleotide sequences use the following command

`python database_curator.py -in (input_file) -n`

if you want to process a nucleotide sequences with optimum length approach use the following command

`python database_curator.py -in (input_file) -n -optimum`

if you want to process a protein sequences use the following command

`python database_curator.py -in (input_file) -p`

if you want to begin your files with a certain name

`python database_curator.py -in (input_file) -p -desired (your_desired_name)`

if you want to process multiple files in the same directory

`python database_curator.py -in (input_file) -p -multi`

if you want to process a fastq file

`python database_curator.py -in (input_file) -n -fastq`


### Any errors please let me know via an e-mail with the subject "database_curator" to eslam.ebrahim@pharma.cu.edu.eg
