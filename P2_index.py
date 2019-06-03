def P2_index(sequ , wwc = False, sst = False, wwy = False, ssy = False, p2 = False ):
        stop_start_c = ['ATG','TAA','TAG','TGA']
        import re
        codon = re.findall('...',str(sequ))
        codon = [i for i in codon if i not in stop_start_c]
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
        P2 = round((WWC + SST) / (WWY + SSY), 3)

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


