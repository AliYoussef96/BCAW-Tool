---
title: 'BCAWT: Automated tool for codon usage bias analysis for molecular evolution'
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

Many measurements have been developed to analyze and study CUB; effective number of codons (ENc) [@Wright1990], codon adaptation index (CAI), relative synonymous codon usage (RSCU) [@Sharp1987] and, translational selection index (P2-index) [@Wang2018]. Also, statistical analysis has been used to investigate the effect of different factors as selection and mutation on shaping CUB such as; Correspondence analysis, Parity Rule 2-plot Analysis and, Neutrality Plot [@Song2017]. BCAWT was developed to analyze such phenomena (Codon Usage Bias) by the aforementioned  measurements.

Various tools are available to analyze and measure CUB, but they lack some important measurements and plots for CUB analysis. What BCAWT does is an automated workflow to study the CUB of any organism's genes by all the measurements and plots mentioned above. Furthermore, using the correlation method to determine the optimal codons described by @Hershberg2009 is implemented for the first time in BCAWT. The tool also includes statistical analysis such as correspondence analysis, correlation analysis, and t-test.
 
# Implementation

BCAWT was developed using python 3.7 with built-in and third-party packages [@Lee2018]. The usage of BCAWT is made to be very easy where users have only to input a fasta format file containing genes to be analyzed, and variety of analyses will be performed.

```python
from BCAWT import BCAWT
BCAWT.BCAW(['Ecoli.fasta'],'save_path',genetic_code_=11,Auto=True)
# processing...
```

The expected outputs from the BCAWT can be divided into three groups. The first one is data in CSV format described in (Table 1), the second group is plots summarized in (Fig 1), and the last one is text files each contains a different result for a different statistical test. The equations used for analyzing CUB in the tool, and the API are reported in the BCAW's [documentation](https://bcaw-tools-documentation.readthedocs.io/en/latest/).

The advantages of BCAWT over existing tools are; 1) the automated workflow, 2) the ability to handle large numbers of genes, 3) the method used to determine optimal codons, named the correlation method, is only available in the BCAWT, 4) visualization and plotting capability, including the creation of violin plots for nucleotide contents, removing the need for other plotting software.

# Output summary

BCAWT returns CSV files containing the CUB indices output (Table 1).

**Table 1: Explanation of the CSV output files from BCAWT.** [(Abbreviations)](https://github.com/AliYoussef96/BCAW-Tool/blob/master/Abbreviations.md) 

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

Furthermore, BCAWT returns 11 plots (Fig 1), enabling an easy interpretation of the results.


**Fig 1: All output plots from BCAWT analysis for coding sequence from Escherichia coli.**

![Fig 1](https://raw.githubusercontent.com/AliYoussef96/BCAW-Tool/master/Plots/All%20plots.jpg)

# References
