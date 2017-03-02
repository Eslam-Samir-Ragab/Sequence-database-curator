# Sequence filteration Program
This program can filter nucleotide and/or protein database from a list of names or sequences (by exact match).
<a href="http://www.freeimagehosting.net/commercial-photography/"><img src="http://i.imgur.com/RL6wMb0.png" alt="commercial photography locations"></a>

## Input:
File containing all the sequences in FASTA format.

## Processing:
It removes specific sequences in your database by your choice.

## Output:
One file of your chosen name.

## Options:
Working on either **protein (p) or nucleotide (n)** databases.

## How to use:
1.	you need to install [python 2.7](https://www.python.org/downloads/) or [python 3](https://www.python.org/downloads/) on your machine.
2. you need to install [Numpy](https://pypi.python.org/pypi/numpy) and [Biopython](http://biopython.org/wiki/Download)
3. you need to install future module by [pip command](https://docs.python.org/3/installing/)
4.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
5.	Open the file “**sequence_filteration.py**” with (python.exe).
  * [Windows](http://stackoverflow.com/a/1527012/7414020)
  * U/Linux : use the command `chmod u+x database_curator.py`
  * Mac : use the command `python sequence_filteration.py`
6.	State your variables and press Enter.

List of options in the program:

| No |    Option   |                                            function                                           |
|:--:|:-----------:|:---------------------------------------------------------------------------------------------:|
|  1 | * -in       | indicate input file, you can use the full path to your file (C:\Users\user\Desktop\dna.fasta) |
|  2 | * -n        | process nucleotide sequences                                                                  |
|  3 | * -p        | process protein sequences                                                                     |
|  4 | * -out      | indicate your output file                                                                     |
|  5 | * -filter   | your list of names or your sequences to be removed                                            |
|  6 | * -flt_mod  | if you want to filter out sequences by names only                                             |
|  7 | -h or -help | show the help for the program                                                                 |

**Required** (*)  -n and -p are mutually exclusive flags


## Examples

if you want to process a nucleotide sequences use the following command

`python sequence_filteration.py -in (input_file) -n -out (output_file) -filter (filter_file) -flt_mode seq`

if you want to process a protein sequences with optimum length approach use the following command

`python sequence_filteration.py -in (input_file) -p -out (output_file) -filter (filter_file) -flt_mode seq`

if you want to process a nucleotide sequences with a file containing list of exact names use the following command

`python sequence_filteration.py -in (input_file) -n -out (output_file) -filter (filter_file) -flt_mode name`

### Any errors please let me know via an e-mail with the subject "database_curator" to eslam.ebrahim@pharma.cu.edu.eg
