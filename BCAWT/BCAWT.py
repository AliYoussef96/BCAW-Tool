def BCAW(main_fasta_file,save_path=str(),ref_fasta_file=None,genetic_code_=1,Auto=False):

    """
    BCAWT ( Bio Codon Analysis Workflow Tool ), it manages a complete workflow to analysis
    the codon usage bias for genes and genomes of any organism..

    Args:

        main_fasta_file (list): list of string of the file's path or file-like object

        save_path (str): absolute path to the directory to save the result in, default = the current directory

        ref_fasta_file (list): list of string of the file's path or file-like object, default = None

        genetic_code_(int): default = 1, The Genetic Codes number described by NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

        Auto (bool): default = False, if ref_fasta_file not None.

    Notes:

        - Auto (bool): should be = True to auto-generate a reference set, when arg (ref_fasta_file) not available ( = None )

    Returns:

       for details see: https://bcaw-tools-documentation.readthedocs.io/en/latest/Table_output.html

    """


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
    import platform
    from pathlib import Path

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    ####Checking if everything is okay
    #####################################
    if ref_fasta_file == None and Auto == False:
        raise ValueError ("No reference genes set specified and Auto is also specified to False")

    if ref_fasta_file != None and Auto == True:
        raise ValueError("Can not specify Auto to be True if reference set genes is specified")

    if ref_fasta_file != None:
        if type(main_fasta_file) != list or type(ref_fasta_file) != list:
            raise TypeError ("main_fasta_file or/and ref_fasta_file have an inappropriate type, should be a list")

    try:
        if ref_fasta_file != None and Auto == False and len(ref_fasta_file) != len(main_fasta_file):
            raise ValueError("The list for both 'main_fasta_file' and 'ref_fasta_file' should be equal")
    except:
        pass

    ############## creating result folder
    ##################################

    if save_path == str():
        directory = Path().absolute()
    else:
        directory = Path(save_path)
    
    print ("Reading Files")


    ##### assign ref_fasta_file if it None
    if ref_fasta_file == None:
        ref_fasta_file = main_fasta_file
    print (main_fasta_file)
    ##########

    ### loop over the input list
    for  input_the_main_fasta_file , input_the_ref_fasta_file in zip (main_fasta_file, ref_fasta_file):


        if Auto == False:
            # append ref. gene set in one list
            print ("Collect the reference genes set")
            all_ref_seq = []
            for seq_ref in SeqIO.parse(input_the_ref_fasta_file, "fasta"):
                lenghth_seq = len(seq_ref.seq) % 3

                if lenghth_seq == 2:
                    seq_main_seq_modifi_ = str(seq_ref.seq[:len(seq_ref.seq) - 2])
                elif lenghth_seq == 1:
                    seq_main_seq_modifi_ = str(seq_ref.seq[:len(seq_ref.seq)-1])
                elif lenghth_seq == 0:
                    seq_main_seq_modifi_ = str(seq_ref.seq)

                all_ref_seq.append( seq_main_seq_modifi_ )

            input_the_ref_fasta_file.seek(0)


    ##############
    ##### first ENc cuz I will used in auto ref
    #############
        print ("ENc")
        print (input_the_main_fasta_file)
        try:
            input_the_main_fasta_file_only = os.path.basename(input_the_main_fasta_file)
            input_the_main_fasta_file_result = input_the_main_fasta_file_only.replace('.fasta', 'ENc.enc')
        except:
            input_the_main_fasta_file_only = os.path.basename(input_the_main_fasta_file.name)
            input_the_main_fasta_file_result = input_the_main_fasta_file_only.replace('.fasta', 'ENc.enc')


        ENc.ENcfilename(input_the_main_fasta_file)

        if type(input_the_main_fasta_file) != str:
            input_the_main_fasta_file.seek(0)   # rest the file object to 0


        ENc_dataframe = pd.read_csv(input_the_main_fasta_file_result, sep='\t')
        ENc_dataframe.drop(['len', 'mo3'], 1, inplace=True)
        ENc_dataframe.rename(index=str, columns={'id': 'gene id', 'eq2Sun': 'ENc'}, inplace=True)


        save_file_name_enc = directory.joinpath(input_the_main_fasta_file_only + "_ENc.csv")

        ENc_dataframe.to_csv(save_file_name_enc, sep=',', index=False, header=True)
        os.remove(input_the_main_fasta_file_result)

        
    #######take the lowest 10% enc as ref.
    #########################################

        if Auto == True:
            print ("Reference gene set")
            list_ref = []

            enc_read_results = pd.read_csv(save_file_name_enc)
            enc_read_results.sort_values(['ENc'], inplace=True)
            data_lowest_enc = int ( round (len(enc_read_results) * 0.10,0) )
            enc_read_results_lowest_enc = enc_read_results.iloc[:data_lowest_enc+1]
            enc_read_results_lowest_enc_gene_id = list(enc_read_results_lowest_enc['gene id'])


            for seq_ref in SeqIO.parse(input_the_main_fasta_file, "fasta"):
                if str(seq_ref.id) in enc_read_results_lowest_enc_gene_id:
                    ### make all seq in ref. gene set divided by 3
                    lenghth_seq = len(str(seq_ref.seq))%3
                    if lenghth_seq == 2:
                        seq_main_seq_modifi_ = str(seq_ref.seq[:len(seq_ref.seq)-2])
                    elif lenghth_seq == 1:
                        seq_main_seq_modifi_ = str(seq_ref.seq[:len(seq_ref.seq)-1])
                    elif lenghth_seq == 0:
                        seq_main_seq_modifi_ = str(seq_ref.seq)
                    list_ref.append(str(seq_main_seq_modifi_))

            if type(input_the_main_fasta_file) != str:
                input_the_main_fasta_file.seek(0)  # rest the file object to 0

            try:
                print (enc_read_results_lowest_enc_gene_id)
            except:
                pass

        print ("______________________________________" + "\n")

    #########################################################


        print("Loading >>> " + str(input_the_main_fasta_file))

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

        for seq_main in SeqIO.parse(input_the_main_fasta_file, "fasta"):

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
            df_ATCG["AROMA"] = GRAVY_AROMO.GRAvy_ARomo(seq_main.seq, genetic_code_, A=True)
            df_ATCG["Gene Length"] = len(str(seq_main.seq))
            df_for_each_file_ATCG = df_for_each_file_ATCG.append(df_ATCG, ignore_index=True, sort=False)

            # if not -auto 1 CAI

            if Auto == False:
                df_CAI = pd.DataFrame()
                df_CAI['gene id'] = [seq_main.id]
                lenghth_seq = len(str(seq_main.seq)) % 3
                if lenghth_seq == 2:
                    seq_main_seq_modifi = str(seq_main.seq[:len(seq_main.seq)-2])
                elif lenghth_seq == 1:
                    seq_main_seq_modifi = str(seq_main.seq[:len(seq_main.seq)-1])
                elif lenghth_seq == 0:
                    seq_main_seq_modifi = str(seq_main.seq)
                df_CAI['CAI'] = [CAI(seq_main_seq_modifi, reference=all_ref_seq,genetic_code = genetic_code_)]
                df_for_each_file_CAI = df_for_each_file_CAI.append(df_CAI, ignore_index=True, sort=False)

            ### cai if == -auto

            elif Auto == True:
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
            ca_rscu = CA_RSCU.CA_RSCU(seq_main.seq,seq_main.id,genetic_code_)
            df_for_each_file_CA_RSCU =pd.concat([ca_rscu, df_for_each_file_CA_RSCU], axis=1)


            print ("Loading >>> " + str(seq_main.id))

        if type(input_the_main_fasta_file) != str:
            input_the_main_fasta_file.seek(0)  # rest the file object to 0


    ############## STEP 2 ( stat )###
    ##################################
    #1 t-test between GC AT
        print ("Loading >>> t-test")
        df_for_each_file_t_test = pd.DataFrame()
        df_for_each_file_t_test.drop(df_for_each_file_t_test.index, inplace=True)  # drop all in df

        t_test_GC_AT = scipy.stats.ttest_rel(df_for_each_file_ATCG['GC'],df_for_each_file_ATCG['AT'],nan_policy = "omit")
        t_test_GC3_AT3 = scipy.stats.ttest_rel(df_for_each_file_ATCG['GC3'], df_for_each_file_ATCG['AT3'], nan_policy="omit")

        try:
            input_the_main_fasta_file_only =os.path.basename(input_the_main_fasta_file)
        except:
            input_the_main_fasta_file_only = os.path.basename(input_the_main_fasta_file.name)


        save_file_name_t_test = directory.joinpath(input_the_main_fasta_file_only +"_t-test.txt")

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
        corr_GC_aromo = scipy.stats.pearsonr(df_for_each_file_ATCG['AROMA'], df_for_each_file_ATCG["GC"])

        #GC3 Vs. GRAvy
        corr_GC3_gravy = scipy.stats.pearsonr(df_for_each_file_ATCG['GRAVY'], df_for_each_file_ATCG["GC3"])
        #GC3 Vs. aromo
        corr_GC3_aromo = scipy.stats.pearsonr(df_for_each_file_ATCG['AROMA'], df_for_each_file_ATCG["GC3"])

        #CAI Vs. gene len
        CAI_Vs_gene_len = scipy.stats.pearsonr(df_for_each_file_ATCG["Gene Length"], df_for_each_file_CAI['CAI'])

        try:
            input_the_main_fasta_file_only =os.path.basename(input_the_main_fasta_file)

        except:
            input_the_main_fasta_file_only = os.path.basename(input_the_main_fasta_file.name)

        save_file_name_Correlation = directory.joinpath(input_the_main_fasta_file_only + "_Correlation.txt")

        with open(save_file_name_Correlation,"w") as f:
            f.write ("Correlation result" + "\n"+
                     "ENc Vs. CAI = " + str (corr_ENc_CAI[0]) + " ,p-value = " + str(corr_ENc_CAI[1])+ "\n"+
                     "ENc Vs. GC = " +  str (corr_ENc_GC[0]) + " ,p-value = " +  str(corr_ENc_GC[1])+ "\n"+
                     "ENc Vs. GC3 = " +  str (corr_ENc_GC3[0]) + " ,p-value = " +  str(corr_ENc_GC3[1])+ "\n"+
                     "CAI Vs. Gene Length = " +  str (corr_cai_gene_len[0]) + " ,p-value = " +  str(corr_cai_gene_len[1])+ "\n"+
                     "GC3 Vs. GC12 = " +  str (corr_GC3_GC12[0]) + " ,p-value = " +  str(corr_GC3_GC12[1])+ "\n"+
                     "GC Vs. GRAVY = "  +  str (corr_GC_gravy[0]) + " ,p-value = " +  str(corr_GC_gravy[1])+ "\n"+
                     "GC Vs. AROMA = "  +  str (corr_GC_aromo[0]) + " ,p-value = " +  str(corr_GC_aromo[1])+ "\n"+
                     "GC3 Vs. GRAVY = "  +  str (corr_GC3_gravy[0]) + " ,p-value = " +  str(corr_GC3_gravy[1])+ "\n"+
                     "GC3 Vs. AROMA = "  +  str (corr_GC3_aromo[0]) + " ,p-value = " +  str(corr_GC3_aromo[1])+ "\n")

    #######Save all DataFrames to csv
    #################################
        print ("Saving Files")
        try:
            input_the_main_fasta_file_only = os.path.basename(input_the_main_fasta_file)
        except:
            input_the_main_fasta_file_only = os.path.basename(input_the_main_fasta_file.name)

        save_file_name_ATCG = directory.joinpath(input_the_main_fasta_file_only + "_ATCG.csv")
        df_for_each_file_ATCG.to_csv(save_file_name_ATCG, sep=',', index=False, header=True)

        save_file_name_CAI = directory.joinpath(input_the_main_fasta_file_only + "_CAI.csv")
        df_for_each_file_CAI.to_csv(save_file_name_CAI, sep=',', index=False, header=True)

        save_file_name_P2 = directory.joinpath(input_the_main_fasta_file_only + "_P2-index.csv")
        df_for_each_file_P2.to_csv(save_file_name_P2, sep=',', index=False, header=True)

        save_file_name_CA_RSCU = directory.joinpath(input_the_main_fasta_file_only + "_CA_RSCU.csv")
        df_for_each_file_CA_RSCU.to_csv(save_file_name_CA_RSCU, sep=',', index=True, header=True)


        ##################
        #Optimal Codon by corr method
        #################
        print ("Loading >>> Optimal Codon")
        Optimal_Codon = Optimal_codon_corr_method.op_corr(save_file_name_enc,save_file_name_CA_RSCU)

        save_file_name_Optimal_Codon = directory.joinpath(input_the_main_fasta_file_only + "_Optimal_Codon.csv")

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
        save_file_name_ENc_CAI_plot = directory.joinpath(input_the_main_fasta_file_only + "_ENc Vs. CAI.png")
        plt.savefig(save_file_name_ENc_CAI_plot)
    #ENc Vs. GC3
        fig2 = plt.figure()
        plt.xlabel("GC3")
        plt.ylabel("ENc")
        plt.title("ENc Vs. GC3")
        plt.scatter(df_for_each_file_ATCG['GC3'],enc_read_results['ENc'],marker ='o',s=10)
        plt.plot(np.unique(df_for_each_file_ATCG['GC3']), np.poly1d(np.polyfit(df_for_each_file_ATCG['GC3'], enc_read_results['ENc'], 1))(np.unique(df_for_each_file_ATCG['GC3'])),'black')

        save_file_name_ENc_GC3_plot = directory.joinpath(input_the_main_fasta_file_only + "_ENc Vs. GC3.png")
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
        save_file_name_Gene_Length_CAI_plot = directory.joinpath(input_the_main_fasta_file_only + "_CAI Vs. Gene Length.png")
        plt.savefig(save_file_name_Gene_Length_CAI_plot)

    #GC Vs. GRAVY
        fig22 = plt.figure()
        plt.xlabel("GRAVY")
        plt.ylabel("GC")
        plt.title("GRAVY Vs. GC")
        plt.scatter(df_for_each_file_ATCG['GRAVY'],df_for_each_file_ATCG['GC'],marker ='o',s=10)
        plt.plot(np.unique(df_for_each_file_ATCG['GRAVY']), np.poly1d(np.polyfit(df_for_each_file_ATCG['GRAVY'], df_for_each_file_ATCG['GC'], 1))(np.unique(df_for_each_file_ATCG['GRAVY'])),'black')

        save_file_name_GC_GRAVY_plot = directory.joinpath(input_the_main_fasta_file_only + "_GRAVY Vs. GC.png")
        plt.savefig(save_file_name_GC_GRAVY_plot)

    #GC Vs. AROMO
        fig33 = plt.figure()
        plt.xlabel("AROMA")
        plt.ylabel("GC")
        plt.title("AROMA Vs. GC")
        plt.scatter(df_for_each_file_ATCG['AROMA'],df_for_each_file_ATCG['GC'],marker ='o',s=10)
        plt.plot(np.unique(df_for_each_file_ATCG['AROMA']), np.poly1d(np.polyfit(df_for_each_file_ATCG['AROMA'], df_for_each_file_ATCG['GC'], 1))(np.unique(df_for_each_file_ATCG['AROMA'])),'black')

        save_file_name_GC_AROMO_plot = directory.joinpath(input_the_main_fasta_file_only + "_AROMA Vs. GC.png")
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

        save_file_name__PR2_plot = directory.joinpath(input_the_main_fasta_file_only + "_PR2-plot.png")
        plt.savefig(save_file_name__PR2_plot)

    #GC12 Vs. GC3
        fig4 = plt.figure()
        plt.xlabel("GC3")
        plt.ylabel("GC12")
        plt.title("GC3 Vs. GC12")
        plt.scatter(df_for_each_file_ATCG['GC3'],df_for_each_file_ATCG['GC12'],marker ='o',s=10)
        plt.plot(np.unique(df_for_each_file_ATCG['GC3']), np.poly1d(np.polyfit(df_for_each_file_ATCG['GC3'], df_for_each_file_ATCG['GC12'], 1))(np.unique(df_for_each_file_ATCG['GC3'])),'black')

        save_file_name_ENc_GC3_plot = directory.joinpath(input_the_main_fasta_file_only + "_GC3 Vs. GC12.png")
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

        save_file_name_ENc_GC_violin = directory.joinpath(input_the_main_fasta_file_only + "_GC violin plot.png")
        plt.savefig(save_file_name_ENc_GC_violin)

    #ACTG3 violin plot
        df_only_ACTG3 = df_for_each_file_ATCG[['A3','T3','C3','G3']]
        fig77 = plt.figure()
        plt.violinplot([df_only_ACTG3['A3'],df_only_ACTG3['T3'],df_only_ACTG3['C3'],df_only_ACTG3['G3']], showmeans=False,showmedians=True)
        plt.title('violin plot for A3, T3, C3 and G3 %')
        plt.ylabel("%")
        #add labels on x-axis
        plt.xticks([1, 2, 3 , 4], ['A3', 'T3', 'C3', 'G3'])

        save_file_name_ACTG3_violin = directory.joinpath(input_the_main_fasta_file_only + "_ACTG3 violin plot.png")
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
        corr_CA_AROMO1 = scipy.stats.pearsonr(read_CA_result_new['axis 1'], df_for_each_file_ATCG['AROMA'])
        #CA_RSCU vs AROMO axis 2
        corr_CA_AROMO2 = scipy.stats.pearsonr(read_CA_result_new['axis 2'], df_for_each_file_ATCG['AROMA'])


        save_file_name_Correlation = directory.joinpath(input_the_main_fasta_file_only + "_Correlation_CA_axis.txt")

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

                     "axis 1 Vs. AROMA = " + str(corr_CA_AROMO1[0]) + " ,p-value = " + str(corr_CA_AROMO1[1]) + "\n" +
                     "axis 2 Vs. AROMA = " + str(corr_CA_AROMO2[0]) + " ,p-value = " + str(corr_CA_AROMO2[1]) + "\n" )

        print ("Results Saved")
