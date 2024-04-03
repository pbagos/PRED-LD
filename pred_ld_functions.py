import dask.dataframe as dd
import numpy as np
import pandas as pd


def Hap_Map_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading Hap Map files ({population}) ...")

    if population not in ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK', 'CHD', 'GIH', "TSI", 'MEX', "ASW"]:
        print("This population is not available in HapMap files. Please select a different population...")
        exit()
    maf_file = f'ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.txt.gz'
    ld_file = f'ref/Hap_Map/ld_chr{chrom}_{population}.txt.gz'
    maf_df = dd.read_csv(maf_file, sep='\s+', blocksize=None)
    # Calculating the Minor Allele Frequency (MAF)
    maf_df['MAF'] = maf_df[['refallele_freq', 'otherallele_freq']].min(axis=1)
    # Filter MAF DataFrame using Dask
    # maf_df = dd.read_csv(maf_file, sep='\s+',blocksize=None, usecols=['rs#', 'chrom', 'pos', 'otherallele_freq'],  )

    # maf_df['het-freq'] = maf_df['het-freq'].astype(float)
    maf_df = maf_df[maf_df['MAF'] >= float(maf_threshold)]

    # Process LD DataFrame using Dask
    ld_df = dd.read_csv(ld_file, blocksize=None, sep='\s+', header=None)

    ld_df = ld_df[(ld_df[3] != ld_df[4]) & (ld_df[6] >= R2_threshold)]

    maf_df = maf_df.rename(columns={'rs#': 'rsID'})

    # Define the new column names
    new_column_names = {
        0: 'pos1',
        1: 'pos2',
        2: 'pop',
        3: 'rsID1',
        4: 'rsID2',
        5: 'Dprime',
        6: 'R2'
    }

    # Rename the columns
    ld_df = ld_df.rename(columns=new_column_names)

    merged_df = dd.merge(ld_df, maf_df, left_on='rsID1', right_on='rsID')
    merged_df = dd.merge(merged_df, maf_df, left_on='rsID2', right_on='rsID')

    merged_df = merged_df[['pos1', 'pos2', 'rsID1', 'rsID2', 'MAF_x', "MAF_y", "R2"]]
    merged_df = merged_df.rename(columns={'MAF_x': 'MAF1', 'MAF_y': 'MAF2'})

    if imp_snp_list:
        final_result = merged_df[merged_df['rsID1'].isin(rs_list) & merged_df['rsID2'].isin(imp_snp_list)]

    else:
        final_result = merged_df[merged_df['rsID1'].isin(rs_list)]
    final_result = final_result.compute()  # Important: This triggers the actual computation
    if final_result.empty:
        print("No SNPs found")
        exit()

    final_result.reset_index(inplace=True, drop=True)
    final_result.to_csv('LD_info_chr' + str(chrom) + '.txt', sep="\t", index=False)
    return final_result


def Hap_Map_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = Hap_Map_LD_info_dask(list(study_df['snp']), chromosome, population, maf_input, r2threshold,
                                      imp_snp_list)

    outputData = pd.merge(outputData, study_df, left_on='rsID1', right_on='snp', how='left')
    outputData['chr'] = chromosome
    # print(outputData.head())
    outputData = outputData.drop(['snp'], axis=1)

    outputData['imputed'] = 0
    #print(outputData.head())
    outputData = outputData.groupby('rsID2').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)

    out_df = pd.DataFrame({
        'snp': outputData['rsID2'],
        'chr': outputData['chr'],
        'pos': outputData['pos2'],
        'beta': outputData['beta'],
        'SE': outputData['SE'],

    })

    pa = 1 - outputData['MAF1'].astype(float)
    pb = 1 - outputData['MAF2'].astype(float)
    pA = outputData['MAF1'].astype(float)
    pB = outputData['MAF2'].astype(float)
    pT = 1 - outputData['MAF1'].astype(float)
    pM = 1 - outputData['MAF2'].astype(float)

    r2_value = outputData['R2']
    D = np.sqrt(r2_value * (pA * pB * pa * pb))
    OR_t = np.exp(outputData['beta'])
    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))
    # f_deriv = (D * (1 - pM)) / (pM * (((1 - pM) * pT - D) * (OR_t - 1) - pM + 1) ** 2)
    # f_deriv = -(D*((2*pM-2)*pT+pM+2*D-1)*OR_t)/((((pM-1)*pT+D)*OR_t+(pM-1)*pT+pM+D-1)*(((pM**2-pM)*pT+D*pM-D)*OR_t+(pM**2-pM)*pT+pM**2+(D-1)*pM+D))
    # Sympy
    # f_deriv = (D*(D - pT*(1 - pM))*(OR_t - 1)/(pM*(-pM + (-D + pT*(1 - pM))*(OR_t - 1) + 1)**2) + D/(pM*(-pM + (-D + pT*(1 - pM))*(OR_t - 1) + 1)))/(D*(OR_t - 1)/(pM*(-pM + (-D + pT*(1 - pM))*(OR_t - 1) + 1)) + 1)
    f_deriv = (-D * (-D + pT * (1 - pM)) * (OR_t - 1) * OR_t / (
                pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1) ** 2) + D * OR_t / (
                           pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1))) / (
                          D * (OR_t - 1) / (pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1)) + 1)

    Var_M = f_deriv ** 2 * var_x ** 2

    out_df['beta'] = np.log(OR_m)
    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['beta'] / out_df['SE']
    out_df['imputed'] = 1
    out_df['R2'] = r2_value
    print(f"Imputed : {sum(out_df['imputed'])} SNPs in chromosome " + str(chromosome))
    # print(out_df.head())
    return out_df


def pheno_Scanner_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    if R2_threshold < 0.8:
        print ("Pheno Scanner data are with a R2 threshold >= 0.8. The R2 threshold will be set to 0.8")
        R2_threshold = 0.8
        
    print("Loading Pheno Scanner files...")

    maf_file = 'ref/pheno_Scanner/1000G.txt'
    ld_file = f'ref/pheno_Scanner/1000G_{population}/1000G_{population}_chr{chrom}.txt.gz'
    population_map = {'EUR': 'eur', 'EAS': 'eas', 'AFR': 'afr', 'AMR': 'amr', 'SAS': 'sas'}

    maf_pop = population_map.get(population, None)
    if maf_pop is None:
        raise ValueError(f"Unsupported population: {population}")

    # Filter MAF DataFrame using Dask
    maf_df = dd.read_csv(maf_file, sep='\s+', blocksize=None, usecols=['hg19_coordinates', 'chr', 'rsid', maf_pop],
                         dtype={maf_pop: 'object'}
                         )
    maf_df = maf_df[(maf_df['chr'] == chrom) & (maf_df[maf_pop] != '-')]
    maf_df[maf_pop] = maf_df[maf_pop].astype(float)
    maf_df = maf_df[1 - maf_df[maf_pop] >= float(maf_threshold)]

    # Process LD DataFrame using Dask
    ld_df = dd.read_csv(ld_file, sep='\s+', blocksize=None, usecols=['ref_hg19_coordinates', 'ref_rsid', 'rsid', 'r2'],
                        dtype={'r2': 'float64'})
    ld_df = ld_df[(ld_df['ref_rsid'] != ld_df['rsid']) & (ld_df['r2'] >= R2_threshold)]

    merged_df = dd.merge(ld_df, maf_df.rename(
        columns={'hg19_coordinates': 'ref_hg19_coordinates', 'rsid': 'ref_rsid', maf_pop: 'MAF1'}), on='ref_rsid')
    merged_df = dd.merge(merged_df, maf_df.rename(columns={maf_pop: 'MAF2'}), on='rsid')

    # Convert to Pandas DataFrame by computing, to finalize and filter based on rs_list
    final_result = merged_df.compute()  # Important: This triggers the actual computation
    # print(final_result.head())
    final_result = final_result.rename(
        columns={'ref_rsid': 'rsID1', 'rsid': 'rsID2', 'ref_hg19_coordinates_x': 'pos1(hg19)',
                 'hg19_coordinates': 'pos2(hg19)', 'r2': 'R2'})
    final_result = final_result[['rsID1', 'pos1(hg19)', 'rsID2', 'pos2(hg19)', 'R2', 'MAF1', 'MAF2']]

    if imp_snp_list:
        final_result = final_result[final_result['rsID2'].isin(rs_list) & final_result['rsID1'].isin(imp_snp_list)]

    else:
        final_result = final_result[final_result['rsID2'].isin(rs_list)]
    if final_result.empty:
        print("No SNPs found")
        exit()
    final_result['pos1(hg19)'] = final_result['pos1(hg19)'].str.split(':').str[1]
    final_result['pos2(hg19)'] = final_result['pos2(hg19)'].str.split(':').str[1]
    final_result.reset_index(inplace=True, drop=True)
    final_result.to_csv('LD_info_chr' + str(chrom) + '.txt', sep="\t", index=False)

    return final_result


def pheno_Scanner_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = pheno_Scanner_LD_info_dask(list(study_df['snp']), chromosome, population, maf_input, r2threshold,
                                            imp_snp_list)
    outputData = pd.merge(outputData, study_df, left_on='rsID2', right_on='snp', how='left')
    outputData['chr'] = chromosome

    outputData = outputData.drop(['snp'], axis=1)
    outputData = outputData.rename(
        columns={"rsID2": "rsID1", "rsID1": "rsID2", "MAF1": "MAF2", "MAF2": "MAF1", "pos1(hg19)": "pos2(hg19)",
                 "pos2(hg19)": "pos1(hg19)"})

    outputData['imputed'] = 0
    # print(outputData.head())
    outputData = outputData.groupby('rsID2').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)

    out_df = pd.DataFrame({
        'snp': outputData['rsID2'],
        'chr': outputData['chr'],
        'pos': outputData['pos2(hg19)'],
        'beta': outputData['beta'],
        'SE': outputData['SE'],

    })

    pa = 1 - outputData['MAF1'].astype(float)
    pb = 1 - outputData['MAF2'].astype(float)
    pA = outputData['MAF1'].astype(float)
    pB = outputData['MAF2'].astype(float)
    pT = outputData['MAF1'].astype(float)
    pM = outputData['MAF2'].astype(float)

    r2_value = outputData['R2']
    D = np.sqrt(r2_value * (pA * pB * pa * pb))
    OR_t = np.exp(outputData['beta'])
    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))
    # f_deriv = (D * (1 - pM)) / (pM * (((1 - pM) * pT - D) * (OR_t - 1) - pM + 1) ** 2)
    # f_deriv = -(D*((2*pM-2)*pT+pM+2*D-1)*OR_t)/((((pM-1)*pT+D)*OR_t+(pM-1)*pT+pM+D-1)*(((pM**2-pM)*pT+D*pM-D)*OR_t+(pM**2-pM)*pT+pM**2+(D-1)*pM+D))
    f_deriv = (-D * (-D + pT * (1 - pM)) * (OR_t - 1) * OR_t / (
                pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1) ** 2) + D * OR_t / (
                           pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1))) / (
                          D * (OR_t - 1) / (pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1)) + 1)

    Var_M = f_deriv ** 2 * var_x ** 2

    out_df['beta'] = np.log(OR_m)
    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['beta'] / out_df['SE']
    out_df['imputed'] = 1
    out_df['R2'] = r2_value

    print(f"Imputed : {sum(out_df['imputed'])} SNPs in chromosome " + str(chromosome))
    # print(out_df.head())
    return out_df


def TOP_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print("Loading TOP-LD files...")
    # Consider converting these to Parquet for better performance
    maf_file = 'ref/TOP_LD/' + population + '/SNV/' + population + '_chr' + str(
        chrom) + '_no_filter_0.2_1000000_info_annotation.csv.gz'
    ld_file = 'ref/TOP_LD/' + population + '/SNV/' + population + '_chr' + str(
        chrom) + '_no_filter_0.2_1000000_LD.csv.gz'

    # Load MAF DataFrame
    maf_df = dd.read_csv(maf_file, blocksize=None, usecols=['Position', 'rsID', 'MAF'])
    # Filter early
    maf_df = maf_df[maf_df['MAF'] >= maf_threshold]

    # Load LD DataFrame
    ld_df = dd.read_csv(ld_file, blocksize=None, usecols=['SNP1', 'SNP2', 'R2'])
    ld_df = ld_df[ld_df['R2'] >= R2_threshold]

    # Merge operations
    # Rename maf_df once and for all
    maf_df = maf_df.rename(columns={'Position': 'SNP', 'rsID': 'rsID', 'MAF': 'MAF'})
    merged_df = dd.merge(ld_df, maf_df.rename(columns={'SNP': 'SNP1', 'rsID': 'rsID1', 'MAF': 'MAF1'}), on='SNP1')
    merged_df = dd.merge(merged_df, maf_df.rename(columns={'SNP': 'SNP2', 'rsID': 'rsID2', 'MAF': 'MAF2'}), on='SNP2')

    # Select and rename desired columns
    final_df = merged_df[['SNP1', 'SNP2', 'R2', 'rsID1', 'rsID2', 'MAF1', 'MAF2']]
    final_df = final_df.rename(columns={'SNP1': 'pos1', 'SNP2': 'pos2'})

    if imp_snp_list:
        result = final_df[final_df['rsID1'].isin(rs_list) & final_df['rsID2'].isin(imp_snp_list)].compute()
    # Compute at the end
    else:
        result = final_df[final_df['rsID1'].isin(rs_list)].compute()

    if result.empty:
        print("No SNPs found")
        exit()

    result.reset_index(inplace=True, drop=True)
    result.to_csv('LD_info_chr' + str(chrom) + '.txt', sep="\t", index=False)
    return result


def TOP_LD_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data

    outputData = TOP_LD_info(list(study_df['snp']), chromosome, population, maf_input, r2threshold, imp_snp_list)

    outputData = pd.merge(outputData, study_df, left_on='rsID1', right_on='snp', how='left')

    outputData = outputData.drop(['snp', 'pos'], axis=1)

    outputData['imputed'] = 0
    outputData = outputData.groupby('rsID2').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)

    out_df = pd.DataFrame({
        'snp': outputData['rsID2'],
        'chr': outputData['chr'],
        'pos': outputData['pos2'],
        'beta': outputData['beta'],
        'SE': outputData['SE'],

    })

    pa = 1 - outputData['MAF1']
    pb = 1 - outputData['MAF2']
    pA = outputData['MAF1']
    pB = outputData['MAF2']
    pT = 1 - outputData['MAF1']
    pM = 1 - outputData['MAF2']

    r2_value = outputData['R2']
    D = np.sqrt(r2_value * (pA * pB * pa * pb))
    OR_t = np.exp(outputData['beta'])
    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))
    # f_deriv = (D * (1 - pM)) / (pM * (((1 - pM) * pT - D) * (OR_t - 1) - pM + 1) ** 2)
    # f_deriv = -(D*((2*pM-2)*pT+pM+2*D-1)*OR_t)/((((pM-1)*pT+D)*OR_t+(pM-1)*pT+pM+D-1)*(((pM**2-pM)*pT+D*pM-D)*OR_t+(pM**2-pM)*pT+pM**2+(D-1)*pM+D))

    # Sympy
    # f_deriv = (D*(D - pT*(1 - pM))*(OR_t - 1)/(pM*(-pM + (-D + pT*(1 - pM))*(OR_t - 1) + 1)**2) + D/(pM*(-pM + (-D + pT*(1 - pM))*(OR_t - 1) + 1)))/(D*(OR_t - 1)/(pM*(-pM + (-D + pT*(1 - pM))*(OR_t - 1) + 1)) + 1)
    f_deriv = (-D * (-D + pT * (1 - pM)) * (OR_t - 1) * OR_t / (
                pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1) ** 2) + D * OR_t / (
                           pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1))) / (
                          D * (OR_t - 1) / (pM * (-pM + (-D + pT * (1 - pM)) * (OR_t - 1) + 1)) + 1)
    Var_M = f_deriv ** 2 * var_x ** 2
    # Var_M = (1/OR_m) ** 2 * var_x ** 2

    out_df['beta'] = np.log(OR_m)
    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['beta'] / out_df['SE']
    out_df['imputed'] = 1
    out_df['R2'] = r2_value

    print(f"Imputed : {sum(out_df['imputed'])} SNPs in chromosome " + str(chromosome))

    return out_df


def process_data(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list):
    final_results_list = []

    study_df = pd.read_csv(file_path, sep="\t")
    chroms = list(set(study_df['chr']))
    ref_panel = ref_file
    # Check if all required columns are present
    required_columns = ['snp', 'chr', 'pos', 'beta', 'SE']  # The required order of columns

    imputed_columns = [col for col in required_columns if col not in required_columns]
    if imputed_columns:
        raise ValueError(
            f" Warning: Check the column names! The columns must be in the following order: snp, chr, pos, beta, SE")

    # Depending on the reference panel...
    if ref_panel == 'TOP_LD':
        for chrom in chroms:
            final_data = TOP_LD_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['beta'] / data['SE']
            data['imputed'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.
            final_df_miss = final_df[final_df['imputed'] == 1]
            final_df_init = final_df[final_df['imputed'] == 0]

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'snp' column.
            final_df_init = final_df_init.drop_duplicates(subset="snp")

            # Concatenate the two DataFrames back together. You might consider resetting the index.
            final_data = pd.concat([final_df_miss, final_df_init]).reset_index(drop=True)
            final_data.to_csv("imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'Pheno_Scanner':
        for chrom in chroms:
            final_data = pheno_Scanner_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['beta'] / data['SE']
            data['imputed'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.
            final_df_miss = final_df[final_df['imputed'] == 1]
            final_df_init = final_df[final_df['imputed'] == 0]

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'snp' column.
            final_df_init = final_df_init.drop_duplicates(subset="snp")

            # Concatenate the two DataFrames back together. You might consider resetting the index.
            final_data = pd.concat([final_df_miss, final_df_init]).reset_index(drop=True)
            final_data.to_csv("imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'Hap_Map':
        for chrom in chroms:
            final_data = Hap_Map_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['beta'] / data['SE']
            data['imputed'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.
            final_df_miss = final_df[final_df['imputed'] == 1]
            final_df_init = final_df[final_df['imputed'] == 0]

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'snp' column.
            final_df_init = final_df_init.drop_duplicates(subset="snp")

            # Concatenate the two DataFrames back together. You might consider resetting the index.
            final_data = pd.concat([final_df_miss, final_df_init]).reset_index(drop=True)
            final_data.to_csv("imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")


    if ref_panel =='all_panels':
        
        print(f"Checking all LD sources")
        for chrom in chroms:
            # For HapMap, we need to take all the panels and merge them...
            if population == 'EUR':
                   
                    pop_hm = "CEU" 
                    final_data_hm  = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                    
                    # pop_hm = "TSI" 
                    # final_data_hm_TSI = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                    # 
                    
                    #final_data_hm = pd.concat([final_data_hm_CEU, final_data_hm_TSI])
                    # Keep the largest R2 value if a snp is common in any of the panels
                    final_data_hm = final_data_hm.groupby('snp').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)
                  
                    
                    
            if population == 'AFR':
                   
                    pop_hm = "YRI" 
                    final_data_hm_YRI  = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                   
                    pop_hm = "MKK" 
                    final_data_hm_MKK = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                   
                    pop_hm = "LWK" 
                    final_data_hm_LWK = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                   
                    pop_hm = "ASW" 
                    final_data_hm_ASW = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
 
                    
                    final_data_hm = pd.concat([final_data_hm_YRI, final_data_hm_MKK,final_data_hm_LWK,final_data_hm_ASW ])
           
            if population == 'EAS':
                   
                    pop_hm = "CHB" 
                    final_data_hm_CHB  = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                    pop_hm = "JPT" 
                    final_data_hm_JPT = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                    pop_hm = "CHD" 
                    final_data_hm_CHD = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)
                    
                    

                    final_data_hm = pd.concat([ final_data_hm_CHB , final_data_hm_JPT, final_data_hm_CHD])
            
            if population == 'SAS':
                  pop_hm = 'GIH'
                  final_data_hm = Hap_Map_process(study_df, r2threshold, pop_hm , maf_input, chrom, imp_snp_list)

            # Keep the largest R2 value if a snp is common in any of the panels
            final_data_hm = final_data_hm.groupby('snp').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)
            final_data_hm['source'] = 'HapMap'
                    
            #final_data_hm = Hap_Map_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_ps = pheno_Scanner_process (study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_tld = TOP_LD_process (study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            final_data_ps['source'] = 'Pheno Scanner'
            final_data_tld['source'] = 'TOP-LD'
           

           # final_data = pd.concat([final_data_hm, final_data_ps, final_data_tld])
            final_data = pd.concat([ final_data_hm, final_data_ps, final_data_tld],ignore_index=True)
       
            # Keep the largest R2 value if a snp is common in any of the panels
            if imp_snp_list == True: 
              final_data = final_data.loc[final_data.groupby('snp')['R2'].idxmax()]
            else : 
              final_data = final_data.groupby('snp').apply(lambda x: x.loc[x['R2'].idxmax()]).reset_index(drop=True)
            print(len(final_data))

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['beta'] / data['SE']
            data['imputed'] = 0
            data['source'] = 'GWAS'
            final_data = pd.concat([final_data, data], ignore_index=True)

            print(f"Total Imputed SNPs: {len(final_data[final_data['imputed'] == 1])} SNPs")

            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            # Separate the DataFrame into two based on the 'imputed' column.
            final_df_miss = final_df[final_df['imputed'] == 1]
            final_df_init = final_df[final_df['imputed'] == 0]

            # Remove duplicates in the 'final_df_init' DataFrame based on the 'snp' column.
            final_df_init = final_df_init.drop_duplicates(subset="snp")

            # Concatenate the two DataFrames back together. You might consider resetting the index.
            final_data = pd.concat([final_df_miss, final_df_init]).reset_index(drop=True)
            final_data.to_csv("imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")
