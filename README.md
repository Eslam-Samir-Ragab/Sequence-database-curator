# Sequence Dereplicator and Database Curator (SDDC) program
This program dereplicates and/or filter nucleotide and/or protein database from a list of names or sequences (by exact match).
<a href="https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/"><img src="https://github.com/Eslam-Samir-Ragab/Sequence-database-curator/blob/master/additionals/description.png"></a>
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

### **The full SDDC commands, Cheat sheet and notes are [here](https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/sddc-program)**

List of options in the program you can download it from [here](https://sites.google.com/pharma.cu.edu.eg/eslam-ibrahim/sddc-program):

<a href="http://www.freeimagehosting.net/commercial-photography/"><img src="http://i.imgur.com/ouwpiBu.png" alt="commercial photography locations"></a>

## Examples

if you want to dereplicate protein sequences use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode derep`

if you want to dereplicate protein sequences with a minimum length = 30 and sequences are in multiple files use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode derep -min_length 30 -multi`

if you want to dereplicate nucleotide sequences with optimum approach and normal protein length = 300 use the following command

`python sddc.py -in (input_file) -n -out (output_file) -mode derep -optimum -prot_length 300`

if you want to filter a protein sequences inclusively by name (i.e. you want to retrieve only seqeunces that you've specified their names) use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode filter -flt_by name -flt_file (filter_file) -approach inclusive`

if you want to filter a protein sequences exclusively by name (i.e. you want to retrieve the seqeunces that aren't present in your filter file) use the following command

`python sddc.py -in (input_file) -p -out (output_file) -mode filter -flt_by name -flt_file (filter_file) -approach exclusive`

if you want to filter a nucleotide sequences by sequence (only exclusive) use the following command

`python sddc.py -in (input_file) -n -out (output_file) -mode filter -flt_by seq -flt_file (filter_file)`

### Any errors please let me know via an e-mail with the subject "database_curator" to eslam.ebrahim@pharma.cu.edu.eg

