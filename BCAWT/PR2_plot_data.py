
def PR2_plot(sequ , o =False, a = False ):
    """Generate data for PR2 plot.

    Args:
        sequ (str): DNA sequence
        o (bool): default = False
        a (bool): default = False

    Returns:
        - ordinate for PR2 plot if arg (o) is True
        - abscissa for PR2 plot if arg (a) is True
    """
    
    import re
    sequ = str(sequ)
    codon = (sequ[i: i + 3] for i in range(0, len(sequ), 3) if len(sequ[i: i + 3]) == 3)
    A3 = 0
    T3 = 0
    G3 = 0
    C3 = 0
    for i_codon in codon:
        if i_codon[2] == 'A':
            A3 += 1
        elif i_codon[2] == 'T':
            T3 += 1
        elif i_codon[2] == 'C':
            C3 += 1
        elif i_codon[2] == 'G':
            G3 += 1

    try:
        ordinate = round ( A3 / (A3 + T3) , 2 )
    except:
        ordinate = 0
    try:
        abscissa = round ( G3 / (G3 + C3) , 2 )
    except:
        abscissa = 0



    if o and o == True:
        return ordinate
    elif a and a == True:
        return abscissa
