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
    
	BCAWT.BCAW('Ecoli','result_folder',fasta=True,Auto=True)
	
    	>> Results Saved
