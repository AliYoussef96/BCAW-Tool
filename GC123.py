#calculating GC1 , GC2 , GC3, GC12
import re
def GC1(sequ):
    from Bio.SeqUtils import GC123
    GC1 = round( GC123(sequ)[1] , 3)
    return GC1

def GC2(sequ):
    from Bio.SeqUtils import GC123
    GC2 = round( GC123(sequ)[2] , 3)
    return GC2

def GC3(sequ):
    from Bio.SeqUtils import GC123
    GC3 = round( GC123(sequ)[3] , 3)
    return GC3

def GC12(sequ):
    from Bio.SeqUtils import GC123
    GC1 = 0
    GC2 = 0
    codon = re.findall('...',str(sequ))
    for i in range(len(codon)):
        if codon[i][0] == 'G' or codon[i][0] == 'C':
            GC1 += 1
        if codon[i][1] == 'G' or codon[i][1] == 'C':
            GC2 += 1
    gc12 = GC1+GC2
    gc12 = round ( gc12 / ( len(codon) * 2 ) , 3)
    return (gc12)


