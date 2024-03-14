# PRED-LD
PRED-LD: A tool for GWAS summary statistics Imputation, using precalculated LD statistics


## Prerequisites
- Pandas  
- NumPy  
- Dask
 

## Installation
- Clone or download this repository to your local machine.
- Ensure that the topld_api executable is correctly placed and has execution permissions. If it's not executable, you can usually change this by running chmod +x topld_api in a terminal.

## Arguments
The script accepts the following command-line arguments:

- --file-path: The path to the input file containing SNP data. The file should be in tab-separated format (TSV,TXT).
- --r2threshold: A float value specifying the R2 threshold for LD filtering.
- --pop: A string indicating the population code to use for LD calculations (EUR, EAS, SAS, AFR).
- --maf: A float value indicating the minor allele frequency (MAF) threshold.
- --ref: A string indicating the LD Reference files (Pheno_Scanner, TOP_LD, Hap_Map)

## Usage
To run the script, navigate to the directory containing the script and execute it with the required arguments. Here is an example command:
```` 
python pred_ld.py --file-path /path/to/your/data.txt --r2threshold 0.8 --pop EUR --maf 0.01 --ref ref_folder_containing_LD_statistics
````

## Example
```` 
python pred_ld.py --file-path TOP_LD_demo.txt --r2threshold 0.8 --pop EUR --maf 0.01 --ref 
````
