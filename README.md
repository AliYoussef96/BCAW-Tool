# BCAW: Automated tool for codon usage bias analysis for molecular evolution

## Statement of Need

There are no tools available enable users to run a whole automated workflow for codon usage bias analysis. Using python 3.7 BCAW Tool ( Bio Codon Analysis Workflow Tool ) was developed to address this problem.
BCAW Tool manages a complete automated workflow to analyze the codon usage bias for genes and genomes of any organism. With minimum coding skills.


For more details about  codon usage bias , and the equations used in BCAWT [see](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Introduction%20to%20codon%20usage%20bias.pdf).


## Dependencies

1- Biopython

2- pandas

3- CAI

4- scipy

5- matplotlib

6- numpy

## Installation Instructions


**Using pip**

```
pip install BCAWT
```

**Note:** Python >=3.7 is required.

## Contribution Guidelines

**Contributions to the software are welcome**

For bugs and suggestions, the most effective way is by raising an issue on the github issue tracker. 
Github allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.

If you wish to contribute some changes to the code then you should submit a [pull request](https://github.com/AliYoussef96/BCAW-Tool/pulls)
How to create a Pull Request? [documentation on pull requests](https://help.github.com/en/articles/about-pull-requests)

## Usage

#### Auto testing

**Note here we try to test the result of BCAW tool and not the modules, for testing the modules in the package use [test.py](https://github.com/AliYoussef96/BCAW-Tool/blob/master/tests/test.py)**

First download fasta file containing the coding sequence ( you can download any fasta file containing gene sequences to be analyzed from [NCBI](https://www.ncbi.nlm.nih.gov/) database).

or just download that file [Test file](https://github.com/AliYoussef96/BCAW-Tool/blob/master/tests/Ecoli.fasta)

then run ( It will automatically run a test on the results files ):

```
from BCAWT import BCAWT_auto_test
BCAWT_auto_test.auto_test(["Ecoli.fasta"])
BCAWT_auto_test.auto_check_files()
>> test is completed 'successfully'
```

#### Main Usage

```
from BCAWT import BCAWT
BCAWT.BCAW(['Ecoli.fasta'],'result_folder',genetic_code_=11,Auto=True)

```
## Input

```

main_fasta_file (list): list of string of the file's path or file-like object

save_folder_name (str): folder name where the result will be saved

ref_fasta_file (list): list of string of the file's path or file-like object, default = None

Auto (bool): default = False, if ref_fasta_file not None.

genetic_code_ (int) : default = 1, The Genetic Codes number described by [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

```
**Important Note:** BCAW tool expect coding sequences as input and not genes, for more information about what the difference between them you can take a look [here](https://qr.ae/TWt2gE)

#### To obtain such fasta file for a species of interest

Say that the species of interest is Escherichia coli str. K-12 substr. MG1655: 

1- Go to the NCBI database.

2- In the search bar write ( Escherichia coli str. K-12 substr. MG1655, complete genome ).

3- choose one of the results ( depending on what you want in your analysis ).

3- On the right of the page, you will find **send to** option. From **sent to** select **Coding Sequences** then **FASTA nucleotides** Finally, press on **Create File**

For [NCBI Genomes Download (FTP) FAQ](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)

## Output

#### The expected CSV files output

|CSV file name|Description|
|------------|-----------|
| ATCG | contains ; gene id, GC, GC1, GC2, GC3, GC12, AT, AT3    A3, T3, C3, G3, GRAVY, AROMO and, Gene Length |
| CA_RSCU | contains ; each RSCU result for each codon in each genes |
| CA_RSCUcodons | contains ; correspondence analysis first 4 axis for each codon |
| CA_RSCUgenes | contains ; correspondence analysis first 4 axis for each gene |
| CAI | contains ; gene id and CAI index |
| ENc | contains ; gene id and ENc index. |
| P2-index | contains ; gene id and P2 index |
| optimal codons | contains; putative optimal codons detected |


#### All output plots from BCAW tool analysis for coding sequence from Escherichia coli

![Fig 1](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Plots/All%20plots.jpg)


## Documentations

1- An intro to the codon usage bias >> [CUB introduction](https://bcaw-tools-documentation.readthedocs.io/en/latest/intro.html)

2. For more information about the equations used to analyze CUB in the BCAW tool >> [Equations](https://bcaw-tools-documentation.readthedocs.io/en/latest/intro.html#equations-used-for-codon-usage-bias-analysis)

2. For more information about the output plots >> [Plots description](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Plots/Plots%20Description.md).

3. For more information about the abbreviations used >> [Abbreviations table](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Abbreviations.md)
