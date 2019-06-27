# BCAW: Automated tool for codon usage bias analysis for molecular evolution

## Statement of Need

There are no tools available enable users to run a whole automated workflow for codon usage bias analysis. Using python 3.7 BCAW Tool ( Bio Codon Analysis Workflow Tool ) was developed to address this problem.
Expected results from BCAW Tool are; nucleotide content, various codon usage measures (as Effective number of codons, Codon adaptation index, etc.) as well as correlation analysis with graphs. Further, Correspondence Analysis and determination of putative optimal codons are included.
BCAW Tool manages a complete automated workflow to analysis the codon usage bias for genes and genomes of any organism. BCAW Tool is available as executable application, work under Windows operating system, also a source code is available.

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
pip install
```

## Contribution Guidelines
To report bugs or seek support please open an issue on this repository. Contributions to the software are welcome; please open an issue for further discussion.

Usage
=====

**Auto testing**

First download fasta file containing the coding sequence

[Test file](https://github.com/AliYoussef96/BCAWTool/blob/master/Ecoli.fasta)

then run ( It will automatically run a test on the results files ):

```
from BCAWT import BCAWT_auto_test
BCAWT_auto_test.auto_test()
BCAWT_auto_test.auto_check_files()
>> test is completed 'successfully'
```

**Main Usage**

```
from BCAWT import BCAWT
BCAWT.BCAW('Ecoli','result_folder',fasta=True,Auto=True)
>> Results Saved
```
**Input**

input_the_main_fasta_file (str) – fasta file contains DNA sequence 

input_the_ref_fasta_file (str) – fasta file contains reference DNA sequence, default = None

**Output**
![Table 1: Expected Output files](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Table.png)


## Output examples

![GC violin plot: Output examples](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Escherichia%20coli%20str.%20K-12%20substr.%20MG1655.fasta_GC%20violin%20plot.png)

GC violin plot: Explain the range of GC, GC1, GC2 and, GC3 content as a normal box plot but also show the probability density of the data at different values for Escherichia coli.

![Correspondence analysis: Output examples](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Escherichia%20coli%20str.%20K-12%20substr.%20MG1655.fasta_CA_RSCU_CA_codos_plot.png)

Correspondence analysis plot for Escherichia coli codons.
