def CA(file):

    """correspondence analysis.

    Args:
        
        file (directory): csv file contains genes' RSCU values

    Returns:
        - csv file contains genes' values for the first 4 axes of the correspondence analysis result
        - csv file contains codons' values for the first 4 axes of the correspondence analysis result
        - plot the genes first 2 axes values of the correspondence analysis result
        - plot the codons first 2 axes values of the correspondence analysis result
    """
    import pandas as pd
    import prince
    import matplotlib.pyplot as plt

    file = str(file)
    df = pd.read_csv(file)
    df.set_index(df.iloc[:,0] , inplace=True)# to make the first column is the index
    df.drop(df.columns[0], axis=1,inplace= True)
    df.replace(0,0.0000001,inplace=True)



    #with prince # make onle CA for 2 axis
    ca = prince.CA(
        n_components=4,
        n_iter=3,
        copy=True,
        check_input=True,
        engine='auto',
        random_state=42
        )
    df.columns.rename('Gene Name', inplace=True)
    df.index.rename('Codons', inplace=True)
    ca = ca.fit(df)

    codons = ca.row_coordinates(df) # for Codons

    genes = ca.column_coordinates(df) #for genes

    #ca.eigenvalues_
    ca.total_inertia_ #total inertia
    ca.explained_inertia_ #inertia for each axis
    inertia = ca.explained_inertia_

    #save information
    file_genes = file.replace(".csv",'')
    file_genes = file_genes + "genes"
    file_genes = file_genes + ".csv"
    genes.rename(columns={genes.columns[0]: 'axis 1', genes.columns[1]: 'axis 2', genes.columns[2]: 'axis 3', genes.columns[3]: 'axis 4'}, inplace=True)
    genes.to_csv(file_genes,sep=',', index=True, header=True) # return csv file for genes ca result



    file_codons = file.replace(".csv",'')
    file_codons = file_codons+ "codons"
    file_codons = file_codons + ".csv"
    codons.rename(columns={codons.columns[0]: 'axis 1', codons.columns[1]: 'axis 2', codons.columns[2]: 'axis 3', codons.columns[3]: 'axis 4'},inplace=True)
    codons.to_csv(file_codons, sep=',', index=True, header=True) # return csv file for codon ca result


    file_inertia = file.replace('.csv','.txt')

    with open(file_inertia, 'a') as f:
        f.write("explained inertia" + "\n")
    for i in range(len(inertia)):
        i_count = i + 1
        with open(file_inertia,'a') as f:
            f.write ("axis " + str(i_count) + " = " + str(inertia[i]) + "\n" )

    with open(file_inertia,'a') as f:
        f.write("Total Inertia = " + str(ca.total_inertia_))

    #plot For genes
    plt.style.use('seaborn-dark-palette')
    fig = plt.figure()

    plt.xlabel("Axis 1")
    plt.ylabel("Axis 2")
    plt.title("CA-plot")
    plt.scatter(genes['axis 1'],genes['axis 2'],s=10,marker ='o')


    plt.axhline(0, color='black', linestyle='-')
    plt.axvline(0, color='black', linestyle='-')


    save_file_name__ca_plot = file + "_CA_gens_plot.png"
    plt.savefig(save_file_name__ca_plot) # return plot file for gene ca result

    #for codons
    plt.style.use('seaborn-dark-palette')
    fig3 = plt.figure()


    plt.xlabel("Axis 1")
    plt.ylabel("Axis 2")
    plt.title("CA-plot")
    plt.scatter(codons['axis 1'],codons['axis 2'], s=10,marker ='o')

    plt.axhline(0, color='black', linestyle='-')
    plt.axvline(0, color='black', linestyle='-')

    if len(codons) < 200:
        for x , y , t in zip(codons['axis 1'],codons['axis 2'] , codons.index.values):
            x = x * (1 + 0.01)
            y = y * (1 + 0.01)
            plt.text(x,y,t)

    file = file.replace('.csv','')
    save_file_name__ca_codons_plot = file + "_CA_codos_plot.png"
    plt.savefig(save_file_name__ca_codons_plot) # return plot file for codon ca result

    read_genes_file = pd.read_csv(file_genes)
    read_genes_file.rename(columns={genes.columns[0]: 'gene id', genes.columns[1]: 'axis 1', genes.columns[2]: 'axis 2'}, inplace=True)



    return read_genes_file


