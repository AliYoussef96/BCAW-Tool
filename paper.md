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

The redundancy in the genetic code means that apart from methionine and tryptophan, an amino acid is encoded by at least two codons. Different codons for the same amino acid are termed synonymous codons. Synonymous codon usage is strongly influenced by evolutionary forces namely, selection and mutation and may vary strongly within or among organisms (Hershberg and Petrov 2008) The preference of specific codons over others contributes to this variation and this phenomenon is called codon usage bias (CUB).

Many measurements have been developed to analysis and study CUB; effective number of codons ( ENc ) (Wright 1990), codon adaptation index ( CAI )  (Sharpl and Li 1987), relative synonymous codon usage ( RSCU ) (Sharp and Li 1987) and, translational selection index ( P2-index ) (Wang et al. 2018). Also, statistical analysis has been used to investigate the effect of different factors as selection and mutation on shaping CUB as; Correspondence analysis (Drosophila 1989), Parity Rule 2 -plot Analysis (Sueoka 2001) and, Neutrality Plot (Song et al. 2017; Wang et al. 2018). BCAW tool was developed to analyze such phenomena ( Codon Usage Bias ) by the aforementioned  measurements.
 
# Implementation

BCAW Tool was developed using python 3.7 with build in and third-party packages (Lee 2018). The usage of BCAW Tool is made to be very easy where users have only to input a fasta format file containing genes to be analyzed, and a bunch of analysis will be performed.

The expected outputs are 23 different files, 8 are in comma-separated values ( CSV ) format, 4 are in text format and, 11 graphs in portable network graphics ( PNG ) format. The equations used for analyzing CUB in the tool and the expected outputs are reported in a tutorial file attached with BCAW tool files.

# References

* Lee. 2018. “Python Implementation of Codon Adaptation Index.” Journal of Open Source Software 3: 905. https://doi.org/10.21105/joss.00905.

* Sharp, Paul M, and Wen-hsiung Li. 1987. “Codon Adaptation Index and Its Potential Applications Nucleic Acids Research” 15 (3): 1281–95.

* Sharpl, Paul M, and Wen-hsiung Li. 1987. “Potential Applications Nucleic Acids Research.” Nucleic Acids Research 15 (3): 1281–95.

* Song, Hui, Jing Liu, Qiuyan Song, Qingping Zhang, Pei Tian, and Zhibiao Nan. 2017. “Comprehensive Analysis of Codon Usage Bias in Seven Epichloë Species and Their Peramine-Coding Genes.” Frontiers in Microbiology 8 (July): 1–12. https://doi.org/10.3389/fmicb.2017.01419.

* Sueoka, Noboru. 2001. “Near Homogeneity of PR2-Bias Fingerprints in the Human Genome and Their Implications in Phylogenetic Analyses.” Molecular Evolution, 469–76. https://doi.org/10.1007/s002390010237.

* Wright, Frank. 1990. “The ‘effective Number of Codons’ Used in a Gene.” Gene 87: 23–29.
