def op_corr(ENc_file_name,RSCU_file_name):
    """
    determine the optimal codons using the correlation method described here: https://doi.org/10.1371/journal.pgen.1000556

    Args:

        ENc_file_name (file): file contains the ENc values for a set of genes
        RSCU_file_name (file): file contains the RSCU values for a set of genes
    Returns:
        DataFrame contains the optimal codons

    """
    #optimal codons by corr. method
    import pandas as pd
    import numpy as np
    import scipy
    from scipy import stats
    from Bio.Data import CodonTable
    standard_table = CodonTable.unambiguous_dna_by_id[1]
    table_codons = standard_table.forward_table

    enc_read_results = pd.read_csv(ENc_file_name)
    enc_read_results.sort_values('gene id', inplace= True)
    enc_read_results.reset_index(drop=True, inplace=True)


    rscu_read_result = pd.read_csv(RSCU_file_name)
    rscu_read_result = rscu_read_result.transpose()
    rscu_read_result.columns = rscu_read_result.iloc[0]
    rscu_read_result = rscu_read_result.iloc[1:]
    rscu_read_result['gene id'] = rscu_read_result.index
    rscu_read_result.sort_values('gene id', inplace=True)
    rscu_read_result.reset_index(drop=True,inplace=True)


    
    df_corr_append = pd.DataFrame()

    for i in rscu_read_result.columns.values:
        if i != 'gene id':
            codon_number_in_family = sum(map((table_codons[i]).__eq__, table_codons.values())) #count values in codon table

            df_corr = pd.DataFrame()
            r = scipy.stats.spearmanr(enc_read_results['ENc'], rscu_read_result[i])
            df_corr["Amino Acid"] = [table_codons[i]]
            df_corr['codons'] = i
            df_corr['r'] = r[0]
            df_corr['p-value'] = r[1]
            df_corr['family_number'] = codon_number_in_family
            df_corr['P_value_/n'] = 0.05 / codon_number_in_family
            #df_corr_append['Optimal Codon'] = ''
            df_corr_append = df_corr_append.append(df_corr, ignore_index=True, sort=False)

    df_corr_append.fillna(0,inplace=True)
    df_corr_append.sort_values(["Amino Acid"], inplace=True)


    df_final_family = pd.DataFrame()
    for i in set(df_corr_append["Amino Acid"]):
        family = df_corr_append['Amino Acid'] == i

        de_family = df_corr_append[family]
        lessthan_p_value = de_family['p-value'] <= de_family['P_value_/n']

        df_lessthan_p_value = de_family[lessthan_p_value]

        min_corr = df_lessthan_p_value['r'] == df_lessthan_p_value['r'].min()

        opetimal_codons = df_lessthan_p_value[min_corr]

        df_final_family = df_final_family.append(opetimal_codons, ignore_index=True, sort=False)

    df_final_family = df_final_family[df_final_family["Amino Acid"] != "W"]
    df_final_family = df_final_family[df_final_family["r"] < 0]
    
    return df_final_family

