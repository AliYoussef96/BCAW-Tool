def GRAvy_ARomo(seq , G = False , A = False):
    """calculating Gravy and Aroma.

    Args:
        seq (str):DNA sequence
        G (bool): default = False
        A (bool): default = False
        

    Returns:
        - Gravy value if arg(G) is True
		
        - Aroma value if arg(A) is True
		
        - None if both args are False

    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.Seq import Seq
    translate_seq = str(seq.translate(table=1))

    protein_seq = translate_seq.replace("*", "")

    protein_seq = ProteinAnalysis(protein_seq)



    AROMO =  protein_seq.aromaticity()
    gravy = protein_seq.gravy()

    if G and G == True:
        return gravy
    elif A and A == True:
        return AROMO

