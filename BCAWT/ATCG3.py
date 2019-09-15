
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
    from itertools import tee

    A3 = 0
    T3 = 0
    G3 = 0
    C3 = 0
    sequ = str(sequ)
    codon, codon_1 = tee(sequ[i: i + 3] for i in range(0, len(sequ), 3) if len(sequ[i: i + 3]) == 3)

    lenght_codon = sum(1 for _ in codon_1)

    for i in codon:
        if A and A == True:
            if i[2] == 'A':
                A3 += 1
        elif T and T == True:
            if i[2] == 'T':
                T3 += 1
        elif C and C == True:
            if i[2] == 'C':
                C3 += 1
        elif G and G == True:
            if i[2] == 'G':
                G3 += 1

    if A and A == True:
        A3 = (A3 / lenght_codon  ) * 100
        return A3
    elif T and T == True:
        T3 = (T3 / lenght_codon ) * 100
        return T3
    elif G and G == True:
        G3 = (G3 / lenght_codon ) * 100
        return G3
    elif C and C == True:
        C3 = (C3 / lenght_codon ) * 100
        return C3


    


