# BCAW: Automated tool for codon usage bias analysis for molecular evolution

## Statement of Need
There are no tools available enable users to run a whole automated workflow for codon usage bias analysis. Using python 3.7 BCAW Tool ( Bio Codon Analysis Workflow Tool ) was developed to address this problem.
Expected results from BCAW Tool are; nucleotide content, various codon usage measures (as Effective number of codons, Codon adaptation index, etc.) as well as correlation analysis with graphs. Further, Correspondence Analysis and determination of putative optimal codons are included.
BCAW Tool manages a complete automated workflow to analysis the codon usage bias for genes and genomes of any organism. BCAW Tool is available as executable application, work under Windows operating system, also a source code is available.
## Dependencies

This software was designed by python 3.7 if it will be used as an executable file on Windows operating system no need for any dependencies. If the source code will be used, these packages must be installed:

1- Biopython

2- pandas

3- CAI

4- scipy

5- matplotlib

6- numpy

## Installation Instructions

1- Download (BCAW Tool.zip) file available at: https://sourceforge.net/projects/bcaw-tool/

2- Extract the files.

3- Within the BCAW Tool folder you expect to find; BCAW tool folder, Source Code folder, Result folder shortcut, BCAW tool.exe shortcut, and pdf tutorial file.

4- Run the tool from BCAW tool.exe shortcut icon.

## Contribution Guidelines
To report bugs or seek support please open an issue on this repository. Contributions to the software are welcome; please open an issue for further discussion.

## How to Use BCAW tool?
Unzip the file you downloaded, open the software from BCAWT shortcut icon then ; 

**The first message from the tool will be:**
```
-s or -m then Enter the name of the file:
```
              
Where user has to input the main ( fasta ) file, here there are two options  to be used;

The first one:
```
-s FileName
```             
FileName refers to a path of ( fasta ) format file contains genes user want to be analyzed. 
Note: Enter filename without ( .fasta )
Example:
```
-s D:\x\y\BCAW Tool\sample_input_files\Escherichia coli str. K-12 substr. MG1655
```
The second one:
```
-m FileName
```
This time FileName refers to a path of ( text ) format file contains paths for ( fasta ) files containing genes to be analyzed for different organisms.

Note: enter filename without ( .txt )

Example:
```
-m D:\x\y\BCAW Tool\sample_input_files\ Genomes
```
**The second message from the tool will be:**
```
-s or -m Enter the name of the reference gene set file or -Auto:
```
Where user has to input the reference genes set ( fasta ) file, here there are three options  to be used; 

The first one:
```
-s FileName_reference
```
FileName_reference, refers to a path of ( fasta ) format file contains the reference genes. 

Note: enter FileName_reference without ( .fasta )

Example:
```
-s D:\x\y\BCAW Tool\sample_input_files\Escherichia coli reference genes
```
The second one:
```
-m FileName reference genes
```
This time FileName_reference refers to a path of ( text ) format file contains paths for ( fasta ) files contains the reference genes.

Note: enter filename without ( .txt )

Example:
```
-m D:\x\y\BCAW Tool\sample_input_files\ Genomes reference genes
 ```
The third one:
```
-Auto
```
Using this option the BCAW tool will build reference genes set and use it.

**The last step is to enter a name where the results will be saved**
```
Folder Name
```
Note: The Folder will be created inside the Result folder shortcut.
Now the user has to wait until the BCAW tool finish the heavy work.

## Expected Output files

![Table 1: Expected Output files](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Table.png)

