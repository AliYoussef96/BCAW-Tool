def GRAvy_ARomo(seq , genetic_code_ = 1, G = False , A = False):
    """calculating Gravy and Aroma for DNA sequence.

    Args:
        seq (str):DNA sequence
        genetic_code_(int): default = 1, The Genetic Codes number described by NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
        G (bool): default = False
        A (bool): default = False
        

    Returns:
        - Gravy value if arg(G) is True

        - Aroma value if arg(A) is True
		
        - None if both args are False

    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.Seq import Seq

    try:
        seq = Seq(seq)
    except:
        pass

    translate_seq = str(seq.translate(table=genetic_code_))

    protein_seq = translate_seq.replace("*", "")

    protein_seq = ProteinAnalysis(protein_seq)



    AROMO =  protein_seq.aromaticity()
    gravy = protein_seq.gravy()

    if G and G == True:
        return gravy
    elif A and A == True:
        return AROMO

