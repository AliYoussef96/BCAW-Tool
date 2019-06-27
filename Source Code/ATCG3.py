def ACTG3(sequ , A = False , T = False , C = False, G = False):
    """
    Calculate A, T, G, and C content at the third position.

    Args:
       
        sequ (str): DNA sequence
        A (bool): default = False
        T (bool): default = False
        C (bool): default = False
        G (bool): default = False

    Returns:
		- A3 content if arg(A) is True		
		- T3 content if arg(T) is True	
		- C3 content if arg(C) is True	
		- G3 content if arg(G) is True	
		- None if all args are False
   
    """
    import re
    A3 = 0
    T3 = 0
    G3 = 0
    C3 = 0
    codon = re.findall('...', str(sequ))
    for i in range(len(codon)):
        if A and A == True:
            if codon[i][2] == 'A':
                A3 += 1
        elif T and T == True:
            if codon[i][2] == 'T':
                T3 += 1
        elif C and C == True:
            if codon[i][2] == 'C':
                C3 += 1
        elif G and G == True:
            if codon[i][2] == 'G':
                G3 += 1

    if A and A == True:
        A3 = (A3 / len(codon)) * 100
        return A3
    elif T and T == True:
        T3 = (T3 / len(codon)) * 100
        return T3
    elif G and G == True:
        G3 = (G3 / len(codon)) * 100
        return G3
    elif C and C == True:
        C3 = (C3 / len(codon)) * 100
        return C3




    


