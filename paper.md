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

Studies have indicated that CUB affects translation elongation and contributes to rapid and accurate protein synthesis (Quax et al. 2015). Furthermore, combined mutation-drift-selection models have been effective in interpreting the implications of CUB (Hershberg and Petrov 2008). 

Many measurements have been developed to analysis and study CUB; effective number of codons ( ENc ) (Wright 1990), codon adaptation index ( CAI )  (Sharpl and Li 1987), relative synonymous codon usage ( RSCU ) (Sharp and Li 1987) and, translational selection index ( P2-index ) (Wang et al. 2018). Also, statistical analysis has been used to investigate the effect of different factors as selection and mutation on shaping CUB as; Correspondence analysis (Drosophila 1989), Parity Rule 2 -plot Analysis (Sueoka 2001) and, Neutrality Plot (Song et al. 2017; Wang et al. 2018)

Various tools are used to analyze and measure CUB, such as CodonW and ACUA, they are widely used but with no graphs outputs and lack some important measurements and plots for CUB analysis that are widely used in previous studies as P2-index (Wang et al. 2018), hydrophilicity and aromaticity (Choudhury, Uddin, and Chakraborty 2018), Parity Rule 2 -plot (Sueoka 2001) and, Neutrality Plot (Chen et al. 2017).

BCAW Tool was developed in such a way to include the indices mentioned above, plots, and additional statistical reports as well as to be automated.
 
# Implementation

BCAW Tool was developed using python 3.7 with build in and third-party packages (Lee 2018).
The usage of BCAW Tool is made to be very easy where users have only to input a fasta format file containing genes to be analyzed, and a bunch of analysis will be performed. For CAI estimation a reference genes set should be provided, two options may be used 1- reference genes set given by the user or 2- an auto option where BCAW Tool generates genes reference set using 10% of genes have the lowest ENc values ( highest biased genes ).

The expected outputs are 23 different files, 8 are in comma-separated values ( CSV ) format, 4 are in text format and, 11 graphs in portable network graphics ( PNG ) format. The equations used in the tool and the expected outputs are reported in a tutorial file attached with BCAW tool files.

# References

* Chen, Ye, Xinxin Li, Xiaojuan Chi, Song Wang, Yanmei Ma, and Jilong Chen. 2017. “Comprehensive Analysis of the Codon Usage Patterns in the Envelope Glycoprotein E2 Gene of the Classical Swine Fever Virus.” PLoS ONE, 1–14. https://doi.org/10.1371/journal. pone.0183646.

* Choudhury, M. N., A. Uddin, and S. Chakraborty. 2018. “Nucleotide Composition and Codon Usage Bias of SRY Gene.” Andrologia 50 (1): 1–11. https://doi.org/10.1111/and.12787.
Drosophila, In. 1989. “Evidence That Mutation Patterns Vary Among Drosophila Transposable Elements.” Journal of Molecular Biology, 843–46.

* Hershberg, Ruth, and Dmitri A. Petrov. 2008. “Selection on Codon Bias.” Annual Review of Genetics 42 (1): 287–99. https://doi.org/10.1146/annurev.genet.42.110807.091442.

* Lee. 2018. “Python Implementation of Codon Adaptation Index.” Journal of Open Source Software 3: 905. https://doi.org/10.21105/joss.00905.

* Quax, Tessa E.F., Nico J. Claassens, Dieter Söll, and John van der Oost. 2015. “Codon Bias as a Means to Fine-Tune Gene Expression.” Molecular Cell 59 (2): 149–61. https://doi.org/10.1016/j.molcel.2015.05.035.

* Sharp, Paul M, and Wen-hsiung Li. 1987. “Codon Adaptation Index and Its Potential Applications Nucleic Acids Research” 15 (3): 1281–95.

* Sharpl, Paul M, and Wen-hsiung Li. 1987. “Potential Applications Nucleic Acids Research.” Nucleic Acids Research 15 (3): 1281–95.

* Song, Hui, Jing Liu, Qiuyan Song, Qingping Zhang, Pei Tian, and Zhibiao Nan. 2017. “Comprehensive Analysis of Codon Usage Bias in Seven Epichloë Species and Their Peramine-Coding Genes.” Frontiers in Microbiology 8 (July): 1–12. https://doi.org/10.3389/fmicb.2017.01419.

* Sueoka, Noboru. 2001. “Near Homogeneity of PR2-Bias Fingerprints in the Human Genome and Their Implications in Phylogenetic Analyses.” Molecular Evolution, 469–76. https://doi.org/10.1007/s002390010237.

* Wang, Liyuan, Huixian Xing, Yanchao Yuan, Xianlin Wang, Muhammad Saeed, Jincai Tao, Wei Feng, Guihua Zhang, Xianliang Song, and Xuezhen Sun. 2018. “Genome-Wide Analysis of Codon Usage Bias in Four Sequenced Cotton Species.” PLoS ONE, 1–17. https://doi.org/10.1371/journal.pone.0194372.

* Wright, Frank. 1990. “The ‘effective Number of Codons’ Used in a Gene.” Gene 87: 23–29.
