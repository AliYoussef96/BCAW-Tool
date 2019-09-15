
def GC1(sequ):
    """Calculates the GC content an the first codon position.

    Args:
        
        sequ (str): DNA sequence

    Returns:
        int: GC1 value
        
    """
    from Bio.SeqUtils import GC123
    GC1 = round( GC123(sequ)[1] , 3)
    return GC1

def GC2(sequ):
    """Calculates the GC content an the 2nd codon position..

    Args:
        
        sequ (str): DNA sequence

    Returns:
        int: GC2 value
        
    """
    from Bio.SeqUtils import GC123
    GC2 = round( GC123(sequ)[2] , 3)
    return GC2

def GC3(sequ):
    """Calculates the GC content an the 3rd codon position..

    Args:
        
        sequ (str): DNA sequence

    Returns:
        int: GC3 value
        
    """
    from Bio.SeqUtils import GC123
    GC3 = round( GC123(sequ)[3] , 3)
    return GC3

def GC12(sequ):
    """Calculates the GC content average at the 1st and 2nd codon positions.

    Args:
        
        sequ (str): DNA sequence

    Returns:
        int: GC12 value
        
    """

    from itertools import tee
    from Bio.SeqUtils import GC123
    GC1 = 0
    GC2 = 0

    sequ = str(sequ)
    codon, codon_1 = tee(sequ[i: i + 3] for i in range(0, len(sequ), 3) if len(sequ[i: i + 3]) == 3)

    lenght_codon = sum(1 for _ in codon_1)

    for i in codon:
        if i[0] == 'G' or i[0] == 'C':
            GC1 += 1
        if i[1] == 'G' or i[1] == 'C':
            GC2 += 1

    gc12 = GC1+GC2

    try:
        gc12 = round ( (gc12 / ( lenght_codon * 2 ) ) , 3)
    except:
        gc12 = 0

    return gc12 * 100


