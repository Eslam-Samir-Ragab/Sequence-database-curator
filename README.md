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
Three files (**curated_seq_only.fasta**, **final.fasta** and **deleted.fasta**) contain only the unique and more complete sequences in your input file depending on your approach.
- *curated_seq_only.fasta* : contains sequences without their annotations.
- *final.fasta* : contains sequences with the exact annotation in your database.
- *deleted.fasta* : contains the names of the deleted seqeunces.

## Options:
1. Working on either **protein (p) or nucleotide (n)** databases.
2. Two approaches (**largest possible length** and **optimum length**).
  * largest possible length approach: gives the longest sequence even if it exceeds the length of your gene.
  * optimum length approach: gives only your gene provided you feed the approximate length of your protein.
![Illustration of both approaches.](https://ppjqaa-ch3302.files.1drv.com/y3p6MyxtnFjVwWErixgUKwFIo5p2TQTrMCdzkWUTBK8yDPWhyqeTJHC8bZwrO1dx1PE9Whj6pKaPSpWg3eiUSNhM59AZBre77KnE7QS95ME1MP7GSne3DjOlJo_0e2JgR_JPLNgR69UHSoZxNPjs0ZY7qEO6utxPfU93PFp7uxMubI/Capture%20%281%29.PNG?psid=1)

## How to use (Python 2.7 or later):
1.	you need to install [python 2.7](https://www.python.org/downloads/) on your machine.
2. you need to install [Numpy](https://pypi.python.org/pypi/numpy) and [Biopython](http://biopython.org/wiki/Download)
3.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
4.	Open the file “**database_curator.2.py**” with (python.exe).
  * [Windows](http://stackoverflow.com/a/1527012/7414020)
  * U/Linux : use the file “**database_curator.linux.2.py**” and the command `chmod u+x database_curator.linux.2.py`
  * Mac : use the file “**database_curator.2.py**” and the command `python database_curator.2.py`
5.	State your variables and press Enter.


## How to use (Python 3):
1.	you need to install [python 3](https://www.python.org/downloads/) on your machine.
2. you need to install [Numpy](https://pypi.python.org/pypi/numpy) and [Biopython](http://biopython.org/wiki/Download)
3.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
4.	Open the file “**database_curator.3.py**” with (python.exe).
  * **[Windows](http://stackoverflow.com/a/1527012/7414020)**
  * **U/Linux** : use the file “**database_curator.linux.3.py**” and the command `chmod u+x database_curator.linux.3.py`
  * **Mac** : use the file “**database_curator.2.py**” and the command `python database_curator.3.py`
5.	State your variables and press Enter.

List of options in the program:

| No |    Option   |                                            function                                           |
|:--:|:-----------:|:---------------------------------------------------------------------------------------------:|
|  1 |     -in     | indicate input file, you can use the full path to your file (C:\Users\user\Desktop\dna.fasta) |
|  2 |      -n     | process nucleotide sequences                                                                  |
|  3 |      -p     | process protein sequences                                                                     |
|  4 |   -desired  | indicate a nomenclature for your output files                                                 |
|  5 |    -multi   | if there are multiple files to process                                                        |
|  6 |   -optimum  | if optimum length approach is wanted'                                                         |
|  7 | -h or -help | show the help for the program                                                                 |

## Examples

if you want to process a nucleotide sequences use the following command

> `python database_curator.2.py -in (input_file) -n`

if you want to process a nucleotide sequences with optimum length approach use the following command

> `python database_curator.2.py -in (input_file) -n -optimum`

if you want to process a protein sequences use the following command

> `python database_curator.2.py -in (input_file) -p`

if you want to begin your files with a certain name

> `python database_curator.2.py -in (input_file) -p -desired (your_desired_name)`

if you want to process multiple files in the same directory

> `python database_curator.2.py -in (input_file) -p -multi`


### Any errors please let me know via an e-mail with the subject "database_curator" to eslam.ebrahim@pharma.cu.edu.eg
