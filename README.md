# Sequence Dereplicator and Database Curator (SDDC) program
This program dereplicates and/or filter nucleotide and/or protein database from a list of names or sequences (by exact match).

## Please, cite: DOI: [10.1007/s00284-017-1327-6](https://link.springer.com/article/10.1007/s00284-017-1327-6)

<a href="https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/github-and-softwares/sddc-program"><img src="https://github.com/Eslam-Samir-Ragab/Sequence-database-curator/blob/master/additionals/description.png"></a>
## How to use:
1.	you need to install [python 2.7](https://www.python.org/downloads/) or [python 3](https://www.python.org/downloads/) on your machine.
2. you need to install [Numpy](https://pypi.python.org/pypi/numpy) and [Biopython](http://biopython.org/wiki/Download)
3. you need to install future module by [pip command](https://docs.python.org/3/installing/)
4.	Click “Clone or download” > “Download ZIP” > extract the downloaded file.
5.	Open the file “**sddc.py**” with (python.exe).
  * [Windows](http://stackoverflow.com/a/1527012/7414020)
  * U/Linux : use the command `chmod u+x sddc.py`
  * Mac : use the command `python sddc.py`
6.	State your variables and press Enter.

### **The full SDDC commands, Cheat sheet and notes are [here](https://github.com/Eslam-Samir-Ragab/Sequence-database-curator/blob/master/additionals/SDDC%20Cheat%20sheet.pdf)**

## *Updates in SDDC v3.0:*
1. Bugs fixes.
2. Usage of -org_order with -kw is updated
3. Exchange FASTA headers mode is now available.

## *Updates in SDDC v2.0:*
1. You can filter the sequences using only keywords (separated by a comma) inclusively or exclusively by adding (-kw) argument to your normal command line.
2. You can get your sequences in their original order after dereplication and/or sequence filtration by adding (-org_order) to your normal command line.

## *Notes:*
* The rate of SDDC as determined using Intel(R) Pentium(R) CPU G630 @ 2.70GHz 2.70 GHz Processor, 4.00 GB RAM, 32-bit Operating System

<a href="https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/github-and-softwares/sddc-program"><img src="https://github.com/Eslam-Samir-Ragab/Sequence-database-curator/blob/master/additionals/rate.png"></a>

* List of options and commands in the program you can download it from [here](https://github.com/Eslam-Samir-Ragab/Sequence-database-curator/blob/master/additionals/SDDC%20Commands.pdf):

<a href="https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/github-and-softwares/sddc-program"><img src="https://github.com/Eslam-Samir-Ragab/Sequence-database-curator/blob/master/additionals/commands.png"></a>

## Examples

if you want to dereplicate protein sequences use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode derep`

if you want to dereplicate protein sequences and preserve the original order of the sequences in the new file use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode derep -org_order`

if you want to dereplicate protein sequences with a minimum length = 30 and sequences are in multiple files use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode derep -min_length 30 -multi`

if you want to dereplicate nucleotide sequences with optimum approach and normal protein length = 300 use the following command

`python sddc.py -in (input_file) -n -out (output_file) -mode derep -optimum -prot_length 300`

if you want to filter a protein sequences inclusively by name (i.e. you want to retrieve only seqeunces that you've specified their names) use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode filter -flt_by name -flt_file (filter_file) -approach inclusive`

if you want to filter a protein sequences inclusively by keyword(s) (i.e. you want to retrieve only seqeunces that you've specified the keywords (separated by a comma) in their names) use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode filter -flt_by name -flt_file (filter_file in csv) -approach inclusive -kw`

if you want to filter a protein sequences exclusively by name (i.e. you want to retrieve the seqeunces that aren't present in your filter file) use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode filter -flt_by name -flt_file (filter_file) -approach exclusive`

if you want to filter a protein sequences exclusively by keyword(s) in their names (i.e. you want to retrieve the seqeunces that certain keywords (separated be a comma) aren't present in your filter file) use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode filter -flt_by name -flt_file (filter_file in csv) -approach exclusive -kw`

if you want to filter a nucleotide sequences by sequence (only exclusive) use the following command

`python sddc.py -in (input_file) -n -out (output_file) -mode filter -flt_by seq -flt_file (filter_file)`

if you want to exchange words in FASTA headers of your protein sequences use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode exchange_headers -ex_file (exchange_file in csv)`

if you want to exchange words in FASTA headers of your nucleotide sequences use the following command

`python sddc.py -in (input_file) -n -out (output_file) -mode exchange_headers -ex_file (exchange_file in csv)`

Example (1)

<a href="https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/github-and-softwares/sddc-program"><img src="https://github.com/Eslam-Samir-Ragab/Sequence-database-curator/blob/master/additionals/example.1.png"></a>

Example (2)

<a href="https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/github-and-softwares/sddc-program"><img src="https://github.com/Eslam-Samir-Ragab/Sequence-database-curator/blob/master/additionals/example.2.png"></a>


### Any errors please send me an email to <eslam.ebrahim@pharma.cu.edu.eg>
## Visit [my website](https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/) for more details, other publications, and contact

