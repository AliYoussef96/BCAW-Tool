Usage
=====

Automated testing
-----------------

First download the `example FASTA file <https://raw.githubusercontent.com/AliYoussef96/BCAW-Tool/master/tests/Ecoli.fasta
>`_ containing a coding sequence then run::

	from BCAWT import BCAWT_auto_test
	
	BCAWT_auto_test.auto_test(["Ecoli.fasta"])
	
	BCAWT_auto_test.auto_check_files()
	
	>> test is completed 'successfully'
	
to automatically run a test on the resulting files.
    
Main Usage
----------------

	from BCAWT import BCAWT
	
	BCAWT.BCAW(['Ecoli.fasta'],'result_folder',genetic_code_=11,Auto=True)

**Important Note:** BCAW Tool expects coding sequences 

Input
------------

- ``main_fasta_file`` (list): list of string of the file's path or file-like object

- ``save_folder_name`` (str): folder name where the result will be saved

- ``ref_fasta_file`` (list): list of string of the file's path or file-like object, default = None

- ``Auto`` (bool): default = False, if ``ref_fasta_file`` is not None.

- ``genetic_code`` (int): default = 1, the genetic code ID from `the NCBI table <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_


To obtain FASTA file for a species of interest
----------------------------------------------

Say that the species of interest is *Escherichia coli* str. K-12 substr. MG1655: 

1. Go to the NCBI's `database <https://ncbi.nlm.nih/gov>`_.

2. In the search bar write ``Escherichia coli str. K-12 substr. MG1655, complete genome``.

3. Choose one of the results (depending on what you want in your analysis).

3. On the write of the page, you will find **send to** option. From **sent to** select **Coding Sequences** then **FASTA nucleotides** Finally, press on **Create File**

For NCBI Genomes Download questions, see their `FAQ <https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/>`_.

