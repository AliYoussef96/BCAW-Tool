Usage
======

**Auto testing**

First download fasta file containing the coding sequence

https://github.com/AliYoussef96/BCAWTool/blob/master/Ecoli.fasta

then run ( It will automatically run a test on the results files ):

	from BCAWT import BCAWT_auto_test
	
	BCAWT_auto_test.auto_test()
    
	BCAWT_auto_test.auto_check_files()
	
	>> test is completed 'successfully'
    
**Main Usage**

	from BCAWT import BCAWT
    
	BCAWT.BCAW('Ecoli','result_folder', genetic_code_ = 11,fasta=True,Auto=True)
	
    	>> Results Saved

**Important Note:** BCAW tool expect coding sequences 

**Input**

input_the_main_fasta_file (str) : fasta file contains DNA sequence 

input_the_ref_fasta_file (str) : fasta file contains reference DNA sequence, default = None. Here Auto is True, to generate automatically reference genes set

genetic_code_ (int) : default = 1, The Genetic Codes number described by [NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

**Important Note:** BCAW tool expect coding sequences

**To obtain such fasta file for a species of interest**

Say that the species of interest is Escherichia coli str. K-12 substr. MG1655: 

1- Go to the NCBI database.

2- In the search bar write ( Escherichia coli str. K-12 substr. MG1655, complete genome ).

3- choose one of the results ( depending on what you want in your analysis ).

3- On the write of the page, you will find **send to** option. From **sent to** select **Coding Sequences** then **FASTA nucleotides** Finally, press on **Create File**

For [NCBI Genomes Download (FTP) FAQ](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)

