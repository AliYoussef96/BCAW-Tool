
#### the expected input are CDS
#input needed
    #1- fasta file for CDS need to be analyze ( -s for single fasta -m for multifily))
    #2- fasta file for refrence set

########main file
def BCAW(input_the_main_fasta_file,save_folder_name,input_the_ref_fasta_file=None,genetic_code_ = 1, fasta = False,txt=False,Auto=False):
    """
    BCAWT ( Bio Codon Analysis Workflow Tool ), it manages a complete workflow to analysis
    the codon usage bias for genes and genomes of any organism..

    Args:

        input_the_main_fasta_file (str): fasta file contains DNA sequence ( don't enter it with .fasta ) or text file contains paths for fasta files ( don't enter it with .txt )

        save_folder_name (str): folder name where the result will be saved

        input_the_ref_fasta_file (str): fasta file contains reference DNA sequence, default = None

        genetic_code_(int): default = 1, The Genetic Codes number described by NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

        fasta (bool): default = False, if the first arg (input_the_main_fasta_file) is text file.

        txt (bool): default = False, if the first arg (input_the_main_fasta_file) is fasta file.

        Auto (bool): default = False, if input_the_ref_fasta_file not None.

    Notes:
        - fasta (bool): should be = True when the first arg (input_the_main_fasta_file) is fasta file.
		
        - txt (bool): should be = True when the first arg (input_the_main_fasta_file) is text file.
		
        - Auto (bool): should be = True to auto-generate a reference set, when arg (input_the_ref_fasta_file) not available ( = None )

    Returns:
       for details see: github.com/AliYoussef96/BCAW-Tool/blob/master/Table.png


    """

    #bad arguments raise errors

    if fasta == True and txt == True:
        raise TypeError("Two file extensions are specified, ( fasta = True, txt = True ). Only one is allowed.")
    if input_the_ref_fasta_file != None and Auto == True:
        raise TypeError("Both genes reference set and Auto are specified  Only one is allowed.")
    if txt == False and fasta == False :
        raise TypeError ("The file extension is not specified, (fasta = False, txt = False )")





    import Bio
    from Bio import SeqIO
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    from Bio.SeqUtils import GC
    from Bio.Data import CodonTable
    from Bio.Alphabet import generic_dna
    import pandas as pd
    from pandas import DataFrame
    from BCAWT import GC123
    from BCAWT import ATCG3
    from BCAWT import ENc
    from CAI import CAI
    import scipy
    from scipy import stats
    from BCAWT import PR2_plot_data
    from BCAWT import P2_index
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    from BCAWT import CA_RSCU
    from BCAWT import CA
    from Bio.Data import CodonTable
    from BCAWT import Optimal_codon_corr_method
    import scipy.optimize as opt
    import warnings
    import time
    from BCAWT import GRAVY_AROMO
    import sys
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    file_name_list = []
    #input_the_main_fasta_file  = input("-s or -m then Enter the name of the file: ")
    if fasta and fasta == True:
        #input_the_main_fasta_file = input_the_main_fasta_file.replace("-s",'')
        input_the_main_fasta_file = input_the_main_fasta_file.lstrip()
        #input_the_main_fasta_file = input_the_main_fasta_file.replace('\\','')
        input_the_main_fasta_file = input_the_main_fasta_file + ".fasta"
        file_name_list.append(input_the_main_fasta_file)
    elif txt and txt == True:
        #input_the_main_fasta_file = input_the_main_fasta_file.replace("-m",'')
        input_the_main_fasta_file = input_the_main_fasta_file.lstrip()
        input_the_main_fasta_file = input_the_main_fasta_file + ".txt"
        with open(input_the_main_fasta_file,'r') as r:
            r = r.readlines()
        for i in r:
            i = i.replace('\n','')
            i = i + ".fasta"
            file_name_list.append(i)


    ############ref file
    file_name_ref_list = []
    if input_the_ref_fasta_file != None and fasta == True and Auto == False:
        input_the_ref_fasta_file = input_the_ref_fasta_file.lstrip()
        input_the_ref_fasta_file = input_the_ref_fasta_file + ".fasta"


        # append ref. gene set in one list
        all_ref_seq = []
        for seq_ref in SeqIO.parse(input_the_ref_fasta_file, "fasta"):
            lenghth_seq = len(str(seq_ref.seq)) % 3
            if lenghth_seq == 2:
                seq_main_seq_modifi_ = str(seq_ref.seq) + 'A'
            elif lenghth_seq == 1:
                seq_main_seq_modifi_ = str(seq_ref.seq) + 'AA'
            elif lenghth_seq == 0:
                seq_main_seq_modifi_ = str(seq_ref.seq)
            all_ref_seq.append( seq_main_seq_modifi_ )


    elif input_the_ref_fasta_file != None and txt == True and Auto == False:
        input_the_ref_fasta_file = input_the_ref_fasta_file.lstrip()
        input_the_ref_fasta_file = input_the_ref_fasta_file + ".txt"
        with open(input_the_ref_fasta_file,'r') as r:
            r = r.readlines()
        for i in r:
            i = i.replace("\n",'')
            i = i + ".fasta"
            file_name_ref_list.append(i)

        # append ref. gene set in one list
        all_ref_seq = []
        for i_fasta_file_name in file_name_ref_list:
            for seq_ref in SeqIO.parse(input_the_ref_fasta_file, "fasta"):
                lenghth_seq = len(str(seq_ref.seq)) % 3
                if lenghth_seq == 2:
                    seq_main_seq_modifi_ = str(seq_ref.seq) + 'A'
                elif lenghth_seq == 1:
                    seq_main_seq_modifi_ = str(seq_ref.seq) + 'AA'
                elif lenghth_seq == 0:
                    seq_main_seq_modifi_ = str(seq_ref.seq)
                all_ref_seq.append(seq_main_seq_modifi_)
    elif Auto and Auto == True and input_the_ref_fasta_file == None:
        list_ref = []

    ############## STEP 1 ( INDEX )###
    ##################################
    #read and open main fasta file.

    while True:
        file_name_exist = 0
        # creat folder where result will be saved
        #directory = save_folder_name.replace('.fasta', '')

        dirname, filename = os.path.split(os.path.realpath(sys.argv[0]))
        directory = os.path.join(dirname, "Result" , save_folder_name )
        directory = directory + "\\"

        if os.path.exists(directory) == False:
            os.makedirs(directory)
            break
        else:
            raise TypeError ("The folder is already exist")

    print ("Reading Files")
    ##############
    ##### first ENc cuz I will used in auto ref
    ####################################################
    # 2 ENc
    for i_file_name in file_name_list:
        print ("ENc")
        print (i_file_name)
        i_file_name_only = os.path.basename(i_file_name)

        ENc.ENcfilename(i_file_name)
        i_file_name_result = i_file_name_only.replace('.fasta', 'ENc.enc')

        ENc_dataframe = pd.read_csv(i_file_name_result, sep='\t')
        ENc_dataframe.drop(['len', 'mo3'], 1, inplace=True)
        ENc_dataframe.rename(index=str, columns={'id': 'gene id', 'eq2Sun': 'ENc'}, inplace=True)

        save_file_name_enc = directory + i_file_name_only + "_ENc.csv"
        ENc_dataframe.to_csv(save_file_name_enc, sep=',', index=False, header=True)
        os.remove(i_file_name_result)
    #######take the lowest 10% enc as ref.
    #########################################
        print ("Reference gene set")
        if Auto and Auto == True:
            enc_read_results = pd.read_csv(save_file_name_enc)
            enc_read_results.sort_values(['ENc'], inplace=True)
            data_lowest_enc = int ( round (len(enc_read_results) * 0.10,0) )
            enc_read_results_lowest_enc = enc_read_results.iloc[:data_lowest_enc+1]
            enc_read_results_lowest_enc_gene_id = list(enc_read_results_lowest_enc['gene id'])
            for seq_ref in SeqIO.parse(i_file_name, "fasta"):
                if str(seq_ref.id) in enc_read_results_lowest_enc_gene_id:
                    ### make all seq in ref. gene set divided by 3

                    lenghth_seq = len(str(seq_ref.seq))%3
                    if lenghth_seq == 2:
                        seq_main_seq_modifi_ = str(seq_ref.seq) + 'A'
                    elif lenghth_seq == 1:
                        seq_main_seq_modifi_ = str(seq_ref.seq) + 'AA'
                    elif lenghth_seq == 0:
                        seq_main_seq_modifi_ = str(seq_ref.seq)
                    list_ref.append(str(seq_main_seq_modifi_))

            try:
                print (enc_read_results_lowest_enc_gene_id)
            except:
                pass

        print ("______________________________________" + "\n")

    #########################################################


    for i_file_name in file_name_list:
        print("Loading >>> " + str(i_file_name))
        print("__________________________________")

        #creat dataframe to collect ACTG information in
        df_for_each_file_ATCG = pd.DataFrame()
        df_for_each_file_ATCG.drop(df_for_each_file_ATCG.index, inplace=True) #drop all in df

        # creat dataframe to collect CAI information in
        df_for_each_file_CAI = pd.DataFrame()
        df_for_each_file_CAI.drop(df_for_each_file_CAI.index, inplace=True)  # drop all in df


        # creat dataframe to collect P2 index information in
        df_for_each_file_P2 = pd.DataFrame()
        df_for_each_file_P2.drop(df_for_each_file_P2.index, inplace=True)  # drop all in df

        # for CA_RSCU
        df_for_each_file_CA_RSCU = pd.DataFrame()
        df_for_each_file_CA_RSCU.drop(df_for_each_file_CA_RSCU.index, inplace=True)  # drop all in df

        # datafreme for PR2- plot
        df_for_each_file_PR2 = pd.DataFrame()
        df_for_each_file_PR2.drop(df_for_each_file_PR2.index, inplace=True)  # drop all in df

        allseq = ''  # for RSCU

        #dataframe for EXp_ENc
        df_exp_enc = pd.DataFrame()
        df_exp_enc.drop(df_exp_enc.index, inplace=True)  # drop all in df


        for seq_main in SeqIO.parse(i_file_name, "fasta"):

            allseq += seq_main.seq  # for RSCU
            allseq = allseq.upper()  # for RSCU

            # 0 ATCG
            df_ATCG = pd.DataFrame()
            df_ATCG['gene id'] = [seq_main.id]
            df_ATCG['GC'] = GC(seq_main.seq)
            df_ATCG['GC1'] = GC123.GC1(seq_main.seq)
            df_ATCG['GC2'] = GC123.GC2(seq_main.seq)
            df_ATCG['GC3'] = GC123.GC3(seq_main.seq)
            df_ATCG['GC12'] = GC123.GC12(seq_main.seq)
            df_ATCG['AT'] = abs(GC(seq_main.seq) - 100)
            df_ATCG['AT3'] = abs(GC123.GC3(seq_main.seq) - 100)
            df_ATCG['A3']  = ATCG3.ACTG3(seq_main.seq , A = True)
            df_ATCG['T3'] = ATCG3.ACTG3(seq_main.seq, T=True)
            df_ATCG['C3'] = ATCG3.ACTG3(seq_main.seq, C=True)
            df_ATCG['G3']  = ATCG3.ACTG3(seq_main.seq , G = True)
            df_ATCG["GRAVY"] = GRAVY_AROMO.GRAvy_ARomo(seq_main.seq, genetic_code_, G=True)
            df_ATCG["AROMO"] = GRAVY_AROMO.GRAvy_ARomo(seq_main.seq, genetic_code_, A=True)
            df_ATCG["Gene Length"] = len(str(seq_main.seq))
            df_for_each_file_ATCG = df_for_each_file_ATCG.append(df_ATCG, ignore_index=True, sort=False)

            # if not -auto 1 CAI

            if Auto == False:
                df_CAI = pd.DataFrame()
                df_CAI['gene id'] = [seq_main.id]
                lenghth_seq = len(str(seq_main.seq)) % 3
                if lenghth_seq == 2:
                    seq_main_seq_modifi = str(seq_main.seq[:-2])
                elif lenghth_seq == 1:
                    seq_main_seq_modifi = str(seq_main.seq[:-1])
                elif lenghth_seq == 0:
                    seq_main_seq_modifi = str(seq_main.seq)
                df_CAI['CAI'] = [CAI(seq_main_seq_modifi, reference=all_ref_seq,genetic_code = genetic_code_)]
                df_for_each_file_CAI = df_for_each_file_CAI.append(df_CAI, ignore_index=True, sort=False)
            ### cai if == -auto
            elif Auto and Auto == True:
                df_CAI = pd.DataFrame()
                df_CAI['gene id'] = [seq_main.id]
                lenghth_seq = len(str(seq_main.seq)) % 3
                if lenghth_seq == 2:
                    seq_main_seq_modifi = str(seq_main.seq[:-2])
                elif lenghth_seq == 1:
                    seq_main_seq_modifi = str(seq_main.seq[:-1])
                elif lenghth_seq == 0:
                    seq_main_seq_modifi = str(seq_main.seq)
                df_CAI['CAI'] = [CAI(seq_main_seq_modifi, reference=list_ref,genetic_code = genetic_code_)]
                df_for_each_file_CAI = df_for_each_file_CAI.append(df_CAI, ignore_index=True, sort=False)

            #4 p2 index
            df_p2 = pd.DataFrame()
            df_p2['gene id'] = [seq_main.id]
            df_p2['WWC'] = [P2_index.P2_index(seq_main.seq, wwc = True)]
            df_p2['SST'] = [P2_index.P2_index(seq_main.seq, sst=True)]
            df_p2['WWY'] = [P2_index.P2_index(seq_main.seq, wwy=True)]
            df_p2['SSY'] = [P2_index.P2_index(seq_main.seq, ssy=True)]
            df_p2['P2'] = [P2_index.P2_index(seq_main.seq, p2=True)]
            df_for_each_file_P2 = df_for_each_file_P2.append(df_p2, ignore_index=True, sort=False)

            # PR2 plot
            df_pr2 = pd.DataFrame()
            df_pr2['gene id'] = [seq_main.id]
            df_pr2['abscissa'] = PR2_plot_data.PR2_plot(seq_main.seq, a = True)
            df_pr2['ordinate'] = PR2_plot_data.PR2_plot(seq_main.seq, o = True)
            df_for_each_file_PR2 = df_for_each_file_PR2.append(df_pr2, ignore_index=True, sort=False)


            #CA_RSCU
            ca_rscu = CA_RSCU.CA_RSCU(seq_main.seq,seq_main.id)
            df_for_each_file_CA_RSCU =pd.concat([ca_rscu, df_for_each_file_CA_RSCU], axis=1)


            print ("Loading >>> " + str(seq_main.id))





    ############## STEP 2 ( stat )###
    ##################################
    #1 t-test between GC AT
        print ("Loading >>> t-test")
        df_for_each_file_t_test = pd.DataFrame()
        df_for_each_file_t_test.drop(df_for_each_file_t_test.index, inplace=True)  # drop all in df

        t_test_GC_AT = scipy.stats.ttest_rel(df_for_each_file_ATCG['GC'],df_for_each_file_ATCG['AT'],nan_policy = "omit")
        t_test_GC3_AT3 = scipy.stats.ttest_rel(df_for_each_file_ATCG['GC3'], df_for_each_file_ATCG['AT3'], nan_policy="omit")

        i_file_name_only =os.path.basename(i_file_name)
        save_file_name_t_test = directory + i_file_name_only +"_t-test.txt"
        with open (save_file_name_t_test,'w') as f:
            f.write ("t-test result" + "\n"+ "GC Vs. AT >>> "+ "t-test = "+ str (t_test_GC_AT[0]) + " ,p-value = " + str (t_test_GC_AT[1]) + "\n"
                     + "GC3 Vs. AT3 >>> " + "t-test = " + str(t_test_GC3_AT3[0]) + " ,p-value = " + str(t_test_GC3_AT3[1]) + "\n")
    ###2 corr
        print("Loading >>> Correlation analysis")
        #sort all dataframes

        enc_read_results = pd.read_csv(save_file_name_enc)
        enc_read_results.sort_values(['gene id'],inplace = True)
        df_for_each_file_CAI.sort_values(['gene id'], inplace=True)
        df_for_each_file_ATCG.sort_values(['gene id'], inplace = True)



        # ENc vs. CAI
        corr_ENc_CAI = scipy.stats.pearsonr(enc_read_results['ENc'],df_for_each_file_CAI['CAI'])
        #ENc vs GC
        corr_ENc_GC = scipy.stats.pearsonr(enc_read_results['ENc'], df_for_each_file_ATCG['GC'])
        #CAI Vs gene length
        corr_cai_gene_len = scipy.stats.pearsonr(df_for_each_file_CAI['CAI'], df_for_each_file_ATCG['Gene Length'])
        #ENc vs. GC3
        corr_ENc_GC3 = scipy.stats.pearsonr(enc_read_results['ENc'], df_for_each_file_ATCG['GC3'])
        #GC12 Vs. GC3
        corr_GC3_GC12 = scipy.stats.pearsonr(df_for_each_file_ATCG['GC3'], df_for_each_file_ATCG['GC12'])

        #GC Vs. GRAvy
        corr_GC_gravy = scipy.stats.pearsonr(df_for_each_file_ATCG['GRAVY'], df_for_each_file_ATCG["GC"])
        #GC Vs. aromo
        corr_GC_aromo = scipy.stats.pearsonr(df_for_each_file_ATCG['AROMO'], df_for_each_file_ATCG["GC"])

        #GC3 Vs. GRAvy
        corr_GC3_gravy = scipy.stats.pearsonr(df_for_each_file_ATCG['GRAVY'], df_for_each_file_ATCG["GC3"])
        #GC3 Vs. aromo
        corr_GC3_aromo = scipy.stats.pearsonr(df_for_each_file_ATCG['AROMO'], df_for_each_file_ATCG["GC3"])

        #CAI Vs. gene len
        CAI_Vs_gene_len = scipy.stats.pearsonr(df_for_each_file_ATCG["Gene Length"], df_for_each_file_CAI['CAI'])

        i_file_name_only =os.path.basename(i_file_name)
        save_file_name_Correlation = directory + i_file_name_only + "_Correlation.txt"
        with open(save_file_name_Correlation,"w") as f:
            f.write ("Correlation result" + "\n"+
                     "ENc Vs. CAI = " + str (corr_ENc_CAI[0]) + " ,p-value = " + str(corr_ENc_CAI[1])+ "\n"+
                     "ENc Vs. GC = " +  str (corr_ENc_GC[0]) + " ,p-value = " +  str(corr_ENc_GC[1])+ "\n"+
                     "ENc Vs. GC3 = " +  str (corr_ENc_GC3[0]) + " ,p-value = " +  str(corr_ENc_GC3[1])+ "\n"+
                     "CAI Vs. Gene Length = " +  str (corr_cai_gene_len[0]) + " ,p-value = " +  str(corr_cai_gene_len[1])+ "\n"+
                     "GC3 Vs. GC12 = " +  str (corr_GC3_GC12[0]) + " ,p-value = " +  str(corr_GC3_GC12[1])+ "\n"+
                     "GC Vs. GRAVY = "  +  str (corr_GC_gravy[0]) + " ,p-value = " +  str(corr_GC_gravy[1])+ "\n"+
                     "GC Vs. AROMO = "  +  str (corr_GC_aromo[0]) + " ,p-value = " +  str(corr_GC_aromo[1])+ "\n"+
                     "GC3 Vs. GRAVY = "  +  str (corr_GC3_gravy[0]) + " ,p-value = " +  str(corr_GC3_gravy[1])+ "\n"+
                     "GC3 Vs. AROMO = "  +  str (corr_GC3_aromo[0]) + " ,p-value = " +  str(corr_GC3_aromo[1])+ "\n")

    #######Save all DataFrames to csv
    #################################
        print ("Saving Files")
        i_file_name_only = os.path.basename(i_file_name)
        save_file_name_ATCG = directory + i_file_name_only + "_ATCG.csv"
        df_for_each_file_ATCG.to_csv(save_file_name_ATCG, sep=',', index=False, header=True)

        save_file_name_CAI = directory + i_file_name_only + "_CAI.csv"
        df_for_each_file_CAI.to_csv(save_file_name_CAI, sep=',', index=False, header=True)

        save_file_name_P2 = directory + i_file_name_only + "_P2-index.csv"
        df_for_each_file_P2.to_csv(save_file_name_P2, sep=',', index=False, header=True)

        #save_file_name_RSCU = directory + i_file_name + "_RSCU.csv"
        #df_for_each_file_RSCU.to_csv(save_file_name_RSCU, sep=',', index=False, header=True)

        save_file_name_CA_RSCU = directory + i_file_name_only + "_CA_RSCU.csv"
        df_for_each_file_CA_RSCU.to_csv(save_file_name_CA_RSCU, sep=',', index=True, header=True)


        ##################
        #Optimal Codon by corr method
        #################
        print ("Loading >>> Optimal Codon")
        Optimal_Codon = Optimal_codon_corr_method.op_corr(save_file_name_enc,save_file_name_CA_RSCU)
        save_file_name_Optimal_Codon = directory + i_file_name_only + "_Optimal_Codon.csv"
        Optimal_Codon.to_csv(save_file_name_Optimal_Codon, sep=',', index=True, header=True)
    #### Step3 Plots  ###############
    #################################
        print (">>> Ploting")
        plt.style.use('seaborn-dark-palette')
    #ENc Vs CAI
        fig = plt.figure()
        plt.xlabel("ENc")
        plt.ylabel("CAI")
        plt.title("ENc Vs. CAI")
        plt.scatter(enc_read_results['ENc'],df_for_each_file_CAI['CAI'],marker ='o',s=10)
        plt.plot(np.unique(enc_read_results['ENc']), np.poly1d(np.polyfit(enc_read_results['ENc'], df_for_each_file_CAI['CAI'], 1))(np.unique(enc_read_results['ENc'])),'black')
        save_file_name_ENc_CAI_plot = directory + i_file_name_only + "_ENc Vs. CAI.png"
        plt.savefig(save_file_name_ENc_CAI_plot)
    #ENc Vs. GC3
        fig2 = plt.figure()
        plt.xlabel("GC3")
        plt.ylabel("ENc")
        plt.title("ENc Vs. GC3")
        plt.scatter(df_for_each_file_ATCG['GC3'],enc_read_results['ENc'],marker ='o',s=10)
        plt.plot(np.unique(df_for_each_file_ATCG['GC3']), np.poly1d(np.polyfit(df_for_each_file_ATCG['GC3'], enc_read_results['ENc'], 1))(np.unique(df_for_each_file_ATCG['GC3'])),'black')

        save_file_name_ENc_GC3_plot = directory + i_file_name_only + "_ENc Vs. GC3.png"
        plt.savefig(save_file_name_ENc_GC3_plot)

    # CAI Vs gene len
        fig55 = plt.figure()
        plt.xlabel("CAI")
        plt.ylabel("Gene Length")
        plt.title("CAI Vs. Gene Length")
        plt.scatter(df_for_each_file_CAI['CAI'], df_for_each_file_ATCG["Gene Length"], marker='o', s=10)
        plt.plot(np.unique(df_for_each_file_CAI['CAI']),
                 np.poly1d(np.polyfit(df_for_each_file_CAI['CAI'], df_for_each_file_ATCG["Gene Length"], 1))(
                     np.unique(df_for_each_file_CAI['CAI'])), 'black')
        save_file_name_Gene_Length_CAI_plot = directory + i_file_name_only + "_CAI Vs. Gene Length.png"
        plt.savefig(save_file_name_Gene_Length_CAI_plot)

    #GC Vs. GRAVY
        fig22 = plt.figure()
        plt.xlabel("GRAVY")
        plt.ylabel("GC")
        plt.title("GRAVY Vs. GC")
        plt.scatter(df_for_each_file_ATCG['GRAVY'],df_for_each_file_ATCG['GC'],marker ='o',s=10)
        plt.plot(np.unique(df_for_each_file_ATCG['GRAVY']), np.poly1d(np.polyfit(df_for_each_file_ATCG['GRAVY'], df_for_each_file_ATCG['GC'], 1))(np.unique(df_for_each_file_ATCG['GRAVY'])),'black')

        save_file_name_GC_GRAVY_plot = directory + i_file_name_only + "_GRAVY Vs. GC.png"
        plt.savefig(save_file_name_GC_GRAVY_plot)

    #GC Vs. AROMO
        fig33 = plt.figure()
        plt.xlabel("AROMO")
        plt.ylabel("GC")
        plt.title("AROMO Vs. GC")
        plt.scatter(df_for_each_file_ATCG['AROMO'],df_for_each_file_ATCG['GC'],marker ='o',s=10)
        plt.plot(np.unique(df_for_each_file_ATCG['AROMO']), np.poly1d(np.polyfit(df_for_each_file_ATCG['AROMO'], df_for_each_file_ATCG['GC'], 1))(np.unique(df_for_each_file_ATCG['AROMO'])),'black')

        save_file_name_GC_AROMO_plot = directory + i_file_name_only + "_AROMO Vs. GC.png"
        plt.savefig(save_file_name_GC_AROMO_plot)
    #PR2-plot
        fig3 = plt.figure()
        plt.xlabel("G3/(G3+C3)")
        plt.ylabel("A3/(A3+T3)")
        plt.title("PR2-plot")
        plt.scatter(df_for_each_file_PR2['abscissa'],df_for_each_file_PR2['ordinate'],marker ='o',s=10)
        plt.xlim(left=0, right=1)
        plt.ylim (bottom=0, top=1)
        plt.vlines(0.5 , ymin=0, ymax=1)
        plt.hlines(0.5, xmin=0, xmax=1)

        save_file_name__PR2_plot = directory + i_file_name_only + "_PR2-plot.png"
        plt.savefig(save_file_name__PR2_plot)

    #GC12 Vs. GC3
        fig4 = plt.figure()
        plt.xlabel("GC3")
        plt.ylabel("GC12")
        plt.title("GC3 Vs. GC12")
        plt.scatter(df_for_each_file_ATCG['GC3'],df_for_each_file_ATCG['GC12'],marker ='o',s=10)
        plt.plot(np.unique(df_for_each_file_ATCG['GC3']), np.poly1d(np.polyfit(df_for_each_file_ATCG['GC3'], df_for_each_file_ATCG['GC12'], 1))(np.unique(df_for_each_file_ATCG['GC3'])),'black')

        save_file_name_ENc_GC3_plot = directory + i_file_name_only + "_GC3 Vs. GC12.png"
        plt.savefig(save_file_name_ENc_GC3_plot)

    #GC0123 violin  plot
        df_only_GC0123 = df_for_each_file_ATCG[['GC','GC1','GC2','GC3']]
        fig66 = plt.figure()
        plt.violinplot([df_only_GC0123['GC'],df_only_GC0123['GC1'],df_only_GC0123['GC2'],df_only_GC0123['GC3']], showmeans=False,showmedians=True)
        plt.title('violin plot for GC%')
        #plt.xlabel("GC%")
        plt.ylabel("%")
        #add labels on x-axis
        plt.xticks([1, 2, 3 , 4], ['GC', 'GC1', 'GC2', 'GC3'])

        save_file_name_ENc_GC_violin = directory + i_file_name_only + "_GC violin plot.png"
        plt.savefig(save_file_name_ENc_GC_violin)

    #ACTG3 violin plot
        df_only_ACTG3 = df_for_each_file_ATCG[['A3','T3','C3','G3']]
        fig77 = plt.figure()
        plt.violinplot([df_only_ACTG3['A3'],df_only_ACTG3['T3'],df_only_ACTG3['C3'],df_only_ACTG3['G3']], showmeans=False,showmedians=True)
        plt.title('violin plot for A3, T3, C3 and G3 %')
        plt.ylabel("%")
        #add labels on x-axis
        plt.xticks([1, 2, 3 , 4], ['A3', 'T3', 'C3', 'G3'])

        save_file_name_ACTG3_violin = directory + i_file_name_only + "_ACTG3 violin plot.png"
        plt.savefig(save_file_name_ACTG3_violin)


        print (">>> CA_RSCU")

    # CA-RSCU
        read_CA_result_new  = CA.CA(save_file_name_CA_RSCU)
        read_CA_result_new.sort_values(by='gene id', inplace=True)

        ###CA_RSCU VS. enc, cai , gc, gc3, gravy and aromo
        ################################################
        ## read_CA_result is returned by CA function

        #CA_RSCU Vs. ENC axis 1
        corr_CA_ENc1 = scipy.stats.pearsonr(read_CA_result_new['axis 1'], enc_read_results['ENc'])
        #CA_RSCU Vs. ENC axis 2
        corr_CA_ENc2 = scipy.stats.pearsonr(read_CA_result_new['axis 2'], enc_read_results['ENc'])

        #CA_RSCU vs CAI axis 1
        corr_CA_CAI1 = scipy.stats.pearsonr(read_CA_result_new['axis 1'], df_for_each_file_CAI['CAI'])
        #CA_RSCU vs CAI axis 2
        corr_CA_CAI2 = scipy.stats.pearsonr(read_CA_result_new['axis 2'], df_for_each_file_CAI['CAI'])

        #CA_RSCU vs GC axis 1
        corr_CA_GC1 = scipy.stats.pearsonr(read_CA_result_new['axis 1'], df_for_each_file_ATCG['GC'])
        #CA_RSCU vs GC axis 2
        corr_CA_GC2 = scipy.stats.pearsonr(read_CA_result_new['axis 2'], df_for_each_file_ATCG['GC'])

        #CA_RSCU vs GC3 axis 1
        corr_CA_GC31 = scipy.stats.pearsonr(read_CA_result_new['axis 1'], df_for_each_file_ATCG['GC3'])
        #CA_RSCU vs GC3 axis 2
        corr_CA_GC32 = scipy.stats.pearsonr(read_CA_result_new['axis 2'], df_for_each_file_ATCG['GC3'])

        #CA_RSCU vs gravy axis 1
        corr_CA_gravy1 = scipy.stats.pearsonr(read_CA_result_new['axis 1'], df_for_each_file_ATCG['GRAVY'])
        #CA_RSCU vs gravy axis 2
        corr_CA_gravy2 = scipy.stats.pearsonr(read_CA_result_new['axis 2'], df_for_each_file_ATCG['GRAVY'])


        #CA_RSCU vs AROMO axis 1
        corr_CA_AROMO1 = scipy.stats.pearsonr(read_CA_result_new['axis 1'], df_for_each_file_ATCG['AROMO'])
        #CA_RSCU vs AROMO axis 2
        corr_CA_AROMO2 = scipy.stats.pearsonr(read_CA_result_new['axis 2'], df_for_each_file_ATCG['AROMO'])


        save_file_name_Correlation = directory + i_file_name_only + "_Correlation_CA_axis.txt"
        with open(save_file_name_Correlation,"w") as f:
            f.write ("Correspondence analysis Correlation result" + "\n"
                     "axis 1 Vs. ENc = " + str (corr_CA_ENc1[0]) + " ,p-value = " + str(corr_CA_ENc1[1])+ "\n" +
                     "axis 2 Vs. ENc = " + str(corr_CA_ENc2[0]) + " ,p-value = " + str(corr_CA_ENc2[1]) + "\n" +

                     "axis 1 Vs. CAI = " + str(corr_CA_CAI1[0]) + " ,p-value = " + str(corr_CA_CAI1[1]) + "\n" +
                     "axis 2 Vs. CAI = " + str(corr_CA_CAI2[0]) + " ,p-value = " + str(corr_CA_CAI2[1]) + "\n" +

                     "axis 1 Vs. GC = " + str(corr_CA_GC1[0]) + " ,p-value = " + str(corr_CA_GC1[1]) + "\n" +
                     "axis 2 Vs. GC = " + str(corr_CA_GC2[0]) + " ,p-value = " + str(corr_CA_GC2[1]) + "\n" +

                     "axis 1 Vs. GC3 = " + str(corr_CA_GC31[0]) + " ,p-value = " + str(corr_CA_GC31[1]) + "\n" +
                     "axis 2 Vs. GC3 = " + str(corr_CA_GC32[0]) + " ,p-value = " + str(corr_CA_GC32[1]) + "\n" +

                     "axis 1 Vs. GRAVY = " + str(corr_CA_gravy1[0]) + " ,p-value = " + str(corr_CA_gravy1[1]) + "\n" +
                     "axis 2 Vs. GRAVY = " + str(corr_CA_gravy2[0]) + " ,p-value = " + str(corr_CA_gravy2[1]) + "\n" +

                     "axis 1 Vs. AROMO = " + str(corr_CA_AROMO1[0]) + " ,p-value = " + str(corr_CA_AROMO1[1]) + "\n" +
                     "axis 2 Vs. AROMO = " + str(corr_CA_AROMO2[0]) + " ,p-value = " + str(corr_CA_AROMO2[1]) + "\n" )

        print ("Results Saved")



