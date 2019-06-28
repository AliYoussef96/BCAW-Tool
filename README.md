# BCAW: Automated tool for codon usage bias analysis for molecular evolution

## Statement of Need

There are no tools available enable users to run a whole automated workflow for codon usage bias analysis. Using python 3.7 BCAW Tool ( Bio Codon Analysis Workflow Tool ) was developed to address this problem.
Expected results from BCAW Tool are; nucleotide content, various codon usage measures (as Effective number of codons, Codon adaptation index, etc.) as well as correlation analysis with graphs. Further, Correspondence Analysis and determination of putative optimal codons are included.
BCAW Tool manages a complete automated workflow to analysis the codon usage bias for genes and genomes of any organism. BCAW Tool is available as executable application, work under Windows operating system, also a source code is available.

For more details about CUB, and the equations used in BCAWT [see](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Introduction%20to%20codon%20usage%20bias.pdf).


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

## Contribution Guidelines

**Contributions to the software are welcome**

For bugs and suggestions, the most effective way is by raising an issue on the github issue tracker. 
Github allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.

If you wish to contribute some changes to the code then you should submit a [pull request](https://github.com/AliYoussef96/BCAW-Tool/pulls)
How to create a Pull Request? [documentation on pull requests](https://help.github.com/en/articles/about-pull-requests)

## Usage

### Auto testing

First download fasta file containing the coding sequence ( you can download any fasta file containing gene sequences to be analyzed from [NCBI](https://www.ncbi.nlm.nih.gov/) database for example).

[Test file](https://github.com/AliYoussef96/BCAWTool/blob/master/Ecoli.fasta)

then run ( It will automatically run a test on the results files ):

```
from BCAWT import BCAWT_auto_test
BCAWT_auto_test.auto_test()
BCAWT_auto_test.auto_check_files()
>> test is completed 'successfully'
```

### Main Usage

```
from BCAWT import BCAWT
BCAWT.BCAW('Ecoli','result_folder',genetic_code_=11,fasta=True,Auto=True)
>> Results Saved
```
### Input

input_the_main_fasta_file (str) : fasta file contains DNA sequence 

input_the_ref_fasta_file (str) : fasta file contains reference DNA sequence, default = None. Here Auto is True, to generate automatically reference genes set

genetic_code_ (int) : default = 1, The Genetic Codes number described by [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

**Important Note:** BCAW tool expect coding sequences

### To obtain such fasta file for a species of interest

Say that the species of interest is Escherichia coli str. K-12 substr. MG1655: 

1- Go to the NCBI database.

2- In the search bar write ( Escherichia coli str. K-12 substr. MG1655, complete genome ).

3- choose one of the results ( depending on what you want in your analysis ).

3- On the write of the page, you will find **send to** option. From **sent to** select **Coding Sequences** then **FASTA nucleotides** Finally, press on **Create File**

For [NCBI Genomes Download (FTP) FAQ](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)

### Output

![Table 1: Expected Output files](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Table.png)

### Abbreviations
[Abbreviations table](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Abbreviations.md)

### API documentation

[BCAW tool’s documentation](https://bcaw-tools-documentation.readthedocs.io/en/latest/index.html)

## Output examples

![GC violin plot: Output examples](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Escherichia%20coli%20str.%20K-12%20substr.%20MG1655.fasta_GC%20violin%20plot.png)

GC violin plot: Explain the range of GC, GC1, GC2 and, GC3 content as a normal box plot but also show the probability density of the data at different values for Escherichia coli.

![Correspondence analysis: Output examples](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Escherichia%20coli%20str.%20K-12%20substr.%20MG1655.fasta_CA_RSCU_CA_codos_plot.png)

Correspondence analysis plot for Escherichia coli codons.

## Documentations

1- For more information about codon usage bias (CUB) and equations used to analyze CUB in the BCAW tool >> [Introduction to CUB](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Introduction%20to%20codon%20usage%20bias.pdf).

2- For more information about the BCAW tool and the API >> [BCAW tool’s documentation](https://bcaw-tools-documentation.readthedocs.io/en/latest/index.html).

3- For more information about the description output.
