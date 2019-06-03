#A3,C3,G3,T3 for plant genome
#1- cound N3 in each gene / all number of codons in this gene
#2- sum each type of N3 / number of genes ( to get the average of each N3 for all genes )

def ACTG3(sequ , A = False , T = False , C = False, G = False):
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




    


