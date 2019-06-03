# make the ordinate and the abscissa for PR2 plot for plant genomes

def PR2_plot(sequ , o =False, a = False ):
    import re
    sequ = str(sequ)
    codon = re.findall('...',sequ)
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


    ordinate = round ( A3 / (A3 + T3) , 2 )
    abscissa = round ( G3 / (G3 + C3) , 2 )

    if o and o == True:
        return ordinate
    elif a and a == True:
        return abscissa
