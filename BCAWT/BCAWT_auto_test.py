def auto_test(path=str(), test_file=str()):
    '''

    Run a demo test with 23 results ( see the documentation for more details about the expected output )

    Args:

        path: absolute path to the directory to save the result in
        test_file: absolute path to the fasta file that will be tested 

    Returns:
        None

    '''

    if path == str():
        raise TypeError ("Must specify a path to save result in")
    if test_file == str():
        raise TypeError ("Must specify a test_file path to fasta file")

    from BCAWT import BCAWT
    file  =  open(test_file, "r")

    BCAWT.BCAW(main_fasta_file = [file] , save_path= path, Auto=True)




def auto_check_files(path):
    '''
    Check the expected outputs.

    Args:
        path: absolute path to the directory to save the result in

    Returns:

        None
    '''

    import glob  
    
    text = [ i for i in glob.glob(f"{path}/*.txt")]

    csv = [ i for i in glob.glob(f"{path}/*.csv")]

    png = [i for i in glob.glob(f"{path}/*.png")]

    if len(text) == 4 and len(csv) == 8 and len(png) == 11:
        
        print ("Test is completed 'successfully'")
    else:
        print ("Something going wrong please, see the documentation or contact the developer")


