import dask.dataframe as dd
import time
import numpy as np
import os
import argparse
import pandas as pd

start = time.time()


def Hap_Map_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold):
    print("Loading Hap Map files...")
    if population not in ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK', 'CHD', 'GIH', "TSI", 'MEX', "ASW"]:
        print("This population is not available in HapMap files. Please select a different population...")
        exit()
    maf_file = f'ref/Hap_Map/genotype_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.txt'
    ld_file = f'ref/Hap_Map/ld_chr{chrom}_{population}.txt'

    # Filter MAF DataFrame using Dask
    maf_df = dd.read_csv(maf_file, sep='\s+', usecols=['rs#', 'chrom', 'pos', 'het-freq'],
                         dtype={'eur': 'object'}
                         )

    # maf_df['het-freq'] = maf_df['het-freq'].astype(float)
    maf_df = maf_df[maf_df['het-freq'] >= float(maf_threshold)]

    # Process LD DataFrame using Dask
    ld_df = dd.read_csv(ld_file, sep='\s+', header=None)

    ld_df = ld_df[(ld_df[3] != ld_df[4]) & (ld_df[6] >= R2_threshold)]

    maf_df = maf_df.rename(columns={'rs#': 'rsID', 'het-freq': 'MAF'})

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

    final_result = merged_df[merged_df['rsID1'].isin(rs_list)]
    final_result = final_result.compute()  # Important: This triggers the actual computation
    if final_result.empty:
        print("No SNPs found")
        exit()

    final_result.reset_index(inplace=True, drop=True)
    final_result.to_csv('LD_info_chr' + str(chrom) + '.txt', sep="\t", index=False)
    return final_result


def Hap_Map_process(study_df, r2threshold, population, maf_input, chromosome):
    # Fetch LD info data
    outputData = Hap_Map_LD_info_dask(list(study_df['snp']), chromosome, population, maf_input, r2threshold)

    outputData = pd.merge(outputData, study_df, left_on='rsID1', right_on='snp', how='left')
    outputData['chr'] = chromosome
    print(outputData.head())
    outputData = outputData.drop(['snp'], axis=1)

    outputData['missing'] = 0
    # print(outputData.head())
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
    pT = outputData['MAF1'].astype(float)
    pM = outputData['MAF2'].astype(float)

    r2_value = outputData['R2']
    D = np.sqrt(r2_value * (pA * pB * pa * pb))
    OR_t = outputData['beta']
    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))
    f_deriv = (D * (1 - pM)) / (pM * (((1 - pM) * pT - D) * (OR_t - 1) - pM + 1) ** 2)
    Var_M = f_deriv ** 2 * var_x ** 2

    out_df['beta'] = OR_m
    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['beta'] / out_df['SE']
    out_df['missing'] = 1
    print(f"Imputed : {sum(out_df['missing'])} SNPs in chromosome " + str(chromosome))
    # print(out_df.head())
    return out_df


def pheno_Scanner_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold):
    if R2_threshold < 0.8:
        raise ValueError("To use Pheno Scanner data, try with an R2 threshold > 0.8")

    print("Loading Pheno Scanner files...")

    maf_file = 'ref/pheno_Scanner/1000G.txt'
    ld_file = f'ref/pheno_Scanner/1000G_{population}_chr{chrom}.txt'
    population_map = {'EUR': 'eur', 'EAS': 'eas', 'AFR': 'afr', 'AMR': 'amr', 'SAS': 'sas'}

    maf_pop = population_map.get(population, None)
    if maf_pop is None:
        raise ValueError(f"Unsupported population: {population}")

    # Filter MAF DataFrame using Dask
    maf_df = dd.read_csv(maf_file, sep='\s+', usecols=['hg19_coordinates', 'chr', 'rsid', maf_pop],
                         dtype={'eur': 'object'}
                         )
    maf_df = maf_df[(maf_df['chr'] == chrom) & (maf_df[maf_pop] != '-')]
    maf_df[maf_pop] = maf_df[maf_pop].astype(float)
    maf_df = maf_df[maf_df[maf_pop] >= float(maf_threshold)]

    # Process LD DataFrame using Dask
    ld_df = dd.read_csv(ld_file, sep='\s+', usecols=['ref_hg19_coordinates', 'ref_rsid', 'rsid', 'r2'])
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

    final_result = final_result[final_result['rsID2'].isin(rs_list)]
    if final_result.empty:
        print("No SNPs found")
        exit()

    final_result.reset_index(inplace=True, drop=True)
    final_result.to_csv('LD_info_chr' + str(chrom) + '.txt', sep="\t", index=False)

    return final_result


# This function now processes large files in chunks using pandas, making it more memory-efficient.

def pheno_Scanner_process(study_df, r2threshold, population, maf_input, chromosome):
    # Fetch LD info data
    outputData = pheno_Scanner_LD_info_dask(list(study_df['snp']), chromosome, population, maf_input, r2threshold)
    outputData = pd.merge(outputData, study_df, left_on='rsID2', right_on='snp', how='left')
    outputData['chr'] = chromosome

    outputData = outputData.drop(['snp'], axis=1)
    outputData = outputData.rename(
        columns={"rsID2": "rsID1", "rsID1": "rsID2", "MAF1": "MAF2", "MAF2": "MAF1", "pos1(hg19)": "pos2(hg19)",
                 "pos2(hg19)": "pos1(hg19)"})

    outputData['missing'] = 0
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
    OR_t = outputData['beta']
    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))
    f_deriv = (D * (1 - pM)) / (pM * (((1 - pM) * pT - D) * (OR_t - 1) - pM + 1) ** 2)
    Var_M = f_deriv ** 2 * var_x ** 2

    out_df['beta'] = OR_m
    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['beta'] / out_df['SE']
    out_df['missing'] = 1
    print(f"Imputed : {sum(out_df['missing'])} SNPs in chromosome " + str(chromosome))
    # print(out_df.head())
    return out_df


def TOP_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold):
    print("Loading TOP-LD files...")
    # Consider converting these to Parquet for better performance
    maf_file = 'ref/TOP_LD/' + population + '_chr' + str(chrom) + '_no_filter_0.2_1000000_info_annotation.csv.gz'
    ld_file = 'ref/TOP_LD/' + population + '_chr' + str(chrom) + '_no_filter_0.2_1000000_LD.csv.gz'

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

    # Compute at the end
    result = final_df[final_df['rsID1'].isin(rs_list)].compute()

    if result.empty:
        print("No SNPs found")
        exit()

    result.reset_index(inplace=True, drop=True)
    result.to_csv('LD_info_chr' + str(chrom) + '.txt', sep="\t", index=False)
    return result


def TOP_LD_process(study_df, r2threshold, population, maf_input, chromosome):
    # Fetch LD info data
    outputData = TOP_LD_info(list(study_df['snp']), chromosome, population, maf_input, r2threshold)

    outputData = pd.merge(outputData, study_df, left_on='rsID1', right_on='snp', how='left')

    outputData = outputData.drop(['snp', 'pos'], axis=1)

    outputData['missing'] = 0
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
    pT = outputData['MAF1']
    pM = outputData['MAF2']

    r2_value = outputData['R2']
    D = np.sqrt(r2_value * (pA * pB * pa * pb))
    OR_t = outputData['beta']
    var_x = outputData['SE']
    OR_m = 1 + ((D * (OR_t - 1)) / (pM * ((1 - pM) + (pT * (1 - pM) - D) * (OR_t - 1))))
    f_deriv = (D * (1 - pM)) / (pM * (((1 - pM) * pT - D) * (OR_t - 1) - pM + 1) ** 2)
    Var_M = f_deriv ** 2 * var_x ** 2

    out_df['beta'] = OR_m
    out_df['SE'] = np.sqrt(Var_M)
    out_df['z'] = out_df['beta'] / out_df['SE']
    out_df['missing'] = 1
    print(f"Imputed : {sum(out_df['missing'])} SNPs in chromosome " + str(chromosome))

    return out_df


def process_data(file_path, r2threshold, population, maf_input, ref_file):
    results = []
    counter = 1
    final_results_list = []
    study_df = pd.read_csv(file_path, sep="\t")
    chroms = list(set(study_df['chr']))
    ref_panel = ref_file
    # Check if all required columns are present
    required_columns = ['snp', 'chr', 'pos', 'beta', 'SE']  # The required order of columns

    missing_columns = [col for col in required_columns if col not in required_columns]
    if missing_columns:
        raise ValueError(
            f" Warning: Check the column names! The columns must be in the following order: snp, chr, pos, beta, SE")

    # Depending on the reference panel...
    if ref_panel == 'TOP_LD':
        for chrom in chroms:
            final_data = TOP_LD_process(study_df, r2threshold, population, maf_input, chrom)

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['beta'] / data['SE']
            data['missing'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            final_df.to_csv("imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'Pheno_Scanner':
        for chrom in chroms:
            final_data = pheno_Scanner_process(study_df, r2threshold, population, maf_input, chrom)

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['beta'] / data['SE']
            data['missing'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            final_df.to_csv("imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")

    if ref_panel == 'Hap_Map':
        for chrom in chroms:
            final_data = Hap_Map_process(study_df, r2threshold, population, maf_input, chrom)

            data = pd.read_csv(file_path, sep="\t")
            data['z'] = data['beta'] / data['SE']
            data['missing'] = 0

            final_data = pd.concat([final_data, data], ignore_index=True)
            print(f"Total : {len(final_data)} SNPs")

            final_data.to_csv("imputation_results_chr" + str(chrom) + ".txt", sep="\t", index=False)

            print("Check 'imputation_results_chr" + str(chrom) + ".txt' for the results")
            print("Check 'LD_info_chr" + str(chrom) + ".txt' for LD information")
            final_results_list.append(final_data)
        if len(chroms) > 1:
            final_df = pd.concat(final_results_list)
            final_df.to_csv("imputation_results_chr_all.txt", sep="\t", index=False)
            print("Check 'imputation_results_chr_all.txt' for results")


def main():
    # Parse arguments
    version = '1.0.0'
    print("---------------------------------------------------------------------------------")
    print("PRED-LD : GWAS Summary Statistics Imputation")
    print("Version " + version + "; April 2024")
    print("Copyright (C) 2024 Pantelis Bagos")
    print("Freely distributed under the GNU General Public Licence (GPLv3)")
    print("---------------------------------------------------------------------------------")

    parser = argparse.ArgumentParser(description="Process data in chunks.")
    parser.add_argument('--file-path', type=str, required=True, help='Input file path')
    parser.add_argument('--r2threshold', type=float, required=True, help='R2 threshold')
    parser.add_argument('--pop', type=str, required=True, help='Population')

    parser.add_argument('--maf', type=float, required=True, help='MAF input value')
    parser.add_argument('--ref', type=str, required=True, help='LD Reference files (pheno_Scanner, TOP_LD, Hap_Map)')

    args = parser.parse_args()

    # Accessing arguments using dashes (convert dashes to underscores)
    file_path = args.file_path
    r2threshold = args.r2threshold
    population = args.pop
    maf_input = args.maf
    ref_file = args.ref
    # chromosome = args.chrom

    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found.")
        return

    process_data(file_path, r2threshold, population, maf_input, ref_file)


if __name__ == "__main__":
    main()

end = time.time()
print(f"Total Time: {end - start} seconds")
