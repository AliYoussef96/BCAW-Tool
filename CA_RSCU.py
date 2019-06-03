def CA_RSCU(allseq,allseq_name):
        from Bio.Data import CodonTable
        from Bio.Seq import Seq
        import re
        from Bio.Alphabet import generic_dna
        import pandas as pd
        from pandas import DataFrame
        stop_start_c = ['ATG','TAA','TAG','TGA']

        stop_start_c = []
        xcodontable = CodonTable.unambiguous_dna_by_id[1]
        ycodontable = xcodontable.forward_table
        zcodontable = []
        qcodontable = []
        result = []
        for i in ycodontable:
            zcodontable.append(ycodontable[i])
            qcodontable.append(i)
        for i in zcodontable:
            if zcodontable.count(i) == 1:
                zcodontable.remove(i)
        RSCU = {}

        allseqstr = str(allseq)
        allseqstr = re.findall('...',allseqstr)
        allseqstr = [ i for i in allseqstr if i not in stop_start_c ]
        qcodontable = [ i for i in qcodontable if i not in stop_start_c]
        dic2 = {}
        allseqstr2 = Seq('', generic_dna)
        for i in allseqstr:
            allseqstr2 += i
        aminoacid2 = allseqstr2.translate(table = 1 , stop_symbol ='')
        aminoacid2 = str(aminoacid2)
        RSCUall = {}

        for ii in allseqstr:
            dic2[ii] = dic2.get(ii,0) + 1
        for k in  qcodontable:
            RSCUall[k] = 0
            if k in dic2:
                try:
                    rscu2 = dic2[k] / ((1/ zcodontable.count(ycodontable[k]))*(aminoacid2.count(ycodontable[k])))
                    RSCUall[k] = round(rscu2,6)

                except ZeroDivisionError:
                    pass

        df = pd.DataFrame(index=pd.Series([i for i in RSCUall]))
        df[allseq_name] = [RSCUall[r] for r in RSCUall ]
        df.drop(['ATG'],0,inplace= True)
        df.sort_index(inplace=True)
        return df
