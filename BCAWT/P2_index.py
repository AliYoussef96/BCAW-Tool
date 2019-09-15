def P2_index(sequ , wwc = False, sst = False, wwy = False, ssy = False, p2 = False ):
    """calculate P2 index.

    Args:
        sequ (str): DNA sequence
        wwc (bool): default = False
        sst (bool): default = False
        wwy (bool): default = False
        ssy (bool): default = False
        p2 (bool): default = False

    Returns:
       - wwc value if arg(wwc) is True
       - sst value if arg(sst) is True
       - wwy value if arg(wwy) is True
       - ssy value if arg(ssy) is True
       - p2 value if arg(p2) is True
    """
    import re
    sequ = str(sequ)
    codon = (sequ[i: i+3] for i in range(0, len(sequ), 3) if len(sequ[i: i+3]) == 3)

    WWC = 0
    SST = 0
    WWY = 0
    SSY = 0
    for i in codon:
        if i == 'AAC' or i == 'TTC' or i == 'ATC' or i == 'TAC':
            WWC += 1
        if i == 'GGT' or i == 'CCT' or i == 'GCT' or i == 'CGT':
            SST += 1
        if i == 'AAT' or i == 'TTT' or i == 'ATT' or i == 'TAT' or i == 'AAC' or i == 'TTC' or i == 'ATC' or i == 'TAC':
            WWY += 1
        if i == 'CCC' or i == 'GGC' or i == 'CGC' or i == 'GCC' or i == 'CCT' or i == 'GGT' or i == 'CGT' or i == 'GCT':
            SSY += 1
    try:
        P2 = round((WWC + SST) / (WWY + SSY), 3)
    except:
        pass

    if wwc and wwc == True:
        return WWC

    elif sst and sst == True:
        return SST

    elif wwy and wwy == True:
        return WWY

    elif ssy and ssy == True:
        return SSY

    elif p2 and p2 == True:
        return P2


