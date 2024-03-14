# PRED-LD
PRED-LD: A tool for GWAS summary statistics Imputation, using precalculated LD statistics
 

## Installation guide
PRED-LD is written in Python (ver. 3.8.2)

1)	Download PRED_LD from: https://github.com/gmanios/PRED-LD

2)	After downloading the .zip folder of MAGE from GitHub, extract it to a working directory. 

3)	Το install the requirements, pip needs to be installed. Download the script for pip, from: https://bootstrap.pypa.io/get-pip.py.

4)	Open a terminal/command prompt, cd to the folder containing the get-pip.py file and run:
    ```
    python get-pip.py
    ```

5)	To install the mentioned requirements with pip, open a terminal/command prompt and run:
    ```
    pip install -r /path/to/requirements.txt
    ```
    
## Arguments
PRED-LD accepts the following command-line arguments:

- --file-path: The path to the input file containing SNP data. The file should be in tab-separated format (TSV,TXT).
- --r2threshold: A float value specifying the R2 threshold for LD filtering.
- --pop: A string indicating the population code to use for LD calculations (EUR, EAS, SAS, AFR).
- --maf: A float value indicating the minor allele frequency (MAF) threshold.
- --ref: A string indicating the LD Reference files (Pheno_Scanner, TOP_LD, Hap_Map)

## Usage
To run PRED-LD, navigate to the directory containing the script and execute it with the required arguments. Here is an example command:
```` 
python pred_ld.py --file-path /path/to/your/data.txt --r2threshold 0.8 --pop EUR --maf 0.01 --ref ref_folder_containing_LD_statistics
````

## Example
```` 
python pred_ld.py --file-path PRED_LD_demo.txt --r2threshold 0.8 --pop EUR --maf 0.01 --ref Pheno_Scanner
````
