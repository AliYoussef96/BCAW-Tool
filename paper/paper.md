---
title: 'BCAW: Automated tool for codon usage bias analysis for molecular evolution'
tags:
  - Codon Usage Analysis
  - Automated Workflow
  - Putative optimal codons
authors:
 - name: Ali Mostafa Anwar
   orcid: 0000-0002-5201-387X
   affiliation: "1"
affiliations:
 - name: Department of Genetics, Faculty of Agriculture, Cairo University, 12613, Cairo, Egypt
   index: 1
date: June 3, 2019
bibliography: paper.bib
---
# Summary

The redundancy in the genetic code means that apart from methionine and tryptophan, an amino acid is encoded by at least two codons. Different codons for the same amino acid are termed synonymous codons. Synonymous codon usage is strongly influenced by evolutionary forces namely, selection and mutation and may vary strongly within or among organisms The preference of specific codons over others contributes to this variation and this phenomenon is called codon usage bias (CUB).

Many measurements have been developed to analysis and study CUB; effective number of codons ( ENc ) [@Wright1990], codon adaptation index ( CAI ), relative synonymous codon usage ( RSCU ) [@Sharp1987] and, translational selection index ( P2-index ) [@Wang2018]. Also, statistical analysis has been used to investigate the effect of different factors as selection and mutation on shaping CUB as; Correspondence analysis, Parity Rule 2-plot Analysis and, Neutrality Plot [@Song2017]. BCAW tool was developed to analyze such phenomena ( Codon Usage Bias ) by the aforementioned  measurements.

Various tools are available to analyze and measure CUB, but they lack some important measurements and plots for CUB analysis. What BCAW tool does is an automated workflow to study the CUB of any organism genes by all the measurements and plots mentioned above. Further, Using correlation method to determine the optimal codons described by [@Hershberg2009] is implemented for the first time in BCAW tool. The tool also includes statistical analysis as correspondence analysis, correlation analysis, and t-test.
 
# Implementation

BCAW Tool was developed using python 3.7 with build in and third-party packages [@Lee2018]. The usage of BCAW Tool is made to be very easy where users have only to input a fasta format file containing genes to be analyzed, and a bunch of analysis will be performed.

```
from BCAWT import BCAWT
BCAWT.BCAW('fasta_file','output_Folder',fasta=True,Auto=True)
```

The expected outputs are 23 different files, 8 are in comma-separated values ( CSV ) format, 4 are in text format and, 11 graphs in portable network graphics ( PNG ) format [(Table)](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Table.png). The equations used for analyzing CUB in the tool and the expected outputs, as well as the API are reported in the BCAWT's [documentation](https://bcaw-tools-documentation.readthedocs.io/en/latest/).

# References
