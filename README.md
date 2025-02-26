# PRED-LD
PRED-LD: A tool for GWAS summary statistics Imputation, using precalculated LD statistics
 
## Web tool
Web tool available at :  [https://compgen.dib.uth.gr/PRED-LD/](https://compgen.dib.uth.gr/PRED-LD/)
  or [195.251.108.198:3839/PRED-LD](http://195.251.108.198:3839/PRED-LD/)

## Demo LD reference panel (Download before running PRED-LD)
[Demo LD ref folder](https://drive.google.com/file/d/1bqz87p0YfxblYrWbFjL4gcq-9fnH8DJP/view?usp=drive_link)


## LD resources 
### HapMap
[HapMap LD Data](https://ftp.ncbi.nlm.nih.gov/hapmap/ld_data/latest/)

[HapMap frequencies](https://ftp.ncbi.nlm.nih.gov/hapmap/frequencies/latest_phaseIII_ncbi_b36/fwd_strand/non-redundant/)

### Pheno Scanner 
[Pheno Scanner database](http://www.phenoscanner.medschl.cam.ac.uk/)

### TOP-LD
[TOP-LD data](http://topld.genetics.unc.edu/downloads/downloads/)


# Population Data Overview
 
## HapMap Populations

The following table lists the populations included in the HapMap project, along with their symbols:
 
| Population | Symbol |
|------------|--------|
| Yoruba in Ibadan, Nigeria | YRI |
| Han Chinese in Beijing, China | CHB |
| Japanese in Tokyo, Japan | JPT |
| CEPH/Utah Collection (NIGMS Human Genetic Cell Repository) | CEU |
| Maasai in Kinyawa, Kenya | MKK |
| Luhya in Webuye, Kenya | LWK |
| Chinese in Metropolitan Denver, CO, USA | CHD |
| Gujarati Indians in Houston, TX, USA | GIH |
| Toscani in Italia | TSI |
| Mexican Ancestry in LA, CA, USA | MXL |
| African Ancestry in SW USA | ASW |
 
 
## Pheno Scanner and TOP-LD Populations

The following table outlines the population symbols as recognized by Pheno Scanner and TOP-LD:
 
| Population | Symbol |
|------------|--------|
| Americans (Only in Pheno Scanner) | AMR |
| South Asians | SAS |
| East Asians | EAS |
| Europeans | EUR |
| Africans | AFR |

 
## Installation guide
PRED-LD is written in Python (ver. 3.8.2)

1)	Clone or download PRED-LD from: https://github.com/pbagos/PRED-LD 
  ```
  git clone  https://github.com/pbagos/PRED-LD
  ```

2)	After downloading the .zip folder of PRED-LD from GitHub, extract it to a working directory. 

3)	Το install the requirements, pip needs to be installed. Download the script for pip, from: https://bootstrap.pypa.io/get-pip.py.

4)	Open a terminal/command prompt, cd to the folder containing the get-pip.py file and run:
    ```
    python get-pip.py
    ```

5)	To install the mentioned requirements with pip, open a terminal/command prompt and run:
    ```
    pip install -r  requirements.txt
    ```
    
## Arguments
PRED-LD accepts the following command-line arguments:

- --file-path: The path to the input file containing SNP data. The file should be in tab-separated format (TSV,TXT) 
- --r2threshold: A float value specifying the R2 threshold for LD filtering 
- --pop: A string indicating the population code to use for LD calculations (EUR, EAS, SAS, AFR, AMR, YRI etc.), depending on the LD reference resource (--ref argument)
- --maf: A float value indicating the minor allele frequency (MAF) threshold
- --ref: A string indicating the LD Reference files (Pheno_Scanner, TOP_LD, Hap_Map, all_panels)
- --imp_list: A filename (.txt) to define specific rsIDs to impute (each SNP has a new line, no header)


## Usage

Your input data should be in a tab-separated text file (TXT format). Ensure the file contains the necessary SNP information and adheres to the specified format:

| snp       | chr | pos       | A1 | A2 | beta       | SE          |
|-----------|-----|-----------|----|----|------------|-------------|
| rs743749  | 22  | 37398195  | A  | G  | -0.6387442 | 9.898344223 |
| rs9306493 | 22  | 45682425  | A  | G  | -0.15022874 | 9.594216875 |
| rs739043  | 22  | 37645230  | G  | A  | -0.05243055 | 9.788226204 |
| rs242885  | 22  | 34423169  | A  | G  | -0.019996628 | 9.449498344 |
| rs5765043 | 22  | 45231883  | G  | A  | -0.07225636 | 9.599864029 |
| rs9625200 | 22  | 27700318  | A  | G  |  0.07320953 | 9.914661823 |
| rs17807317| 22  | 17680519  | C  | A  |  0.5180513 | 9.805693943 |

**Notes:**
- **A1**: Represents the Alternative allele (ALT).
- **A2**: Represents the Reference allele (REF).

To run PRED-LD, navigate to the directory containing the script and execute it with the required arguments. Make sure you have unzipped in the same working directory the ref folder. [Demo LD ref folder](https://drive.google.com/file/d/1mCpiDJZiO9XdBe-6Y0fbXGraF62QqFn5/view?usp=sharing) (Download before running PRED-LD)

## Examples

Here is an example command:
```` 
python pred_ld.py --file-path /path/to/your/data.txt --r2threshold 0.8 --pop EUR --maf 0.01 --ref TOP_LD
````

### Example 1 (Simple Imputation)
```` 
python pred_ld.py --file-path PRED_LD_demo.txt --r2threshold 0.8 --pop EUR --maf 0.01 --ref TOP_LD
````

### Example 2 (Use a list to impute specific rsIDs)
```` 
python pred_ld.py --file-path PRED_LD_demo.txt --r2threshold 0.8 --pop EUR --maf 0.01 --ref TOP_LD --imp_list missing_snps.txt 
````

### Example 3 (Use all panels to perform Imputation)
```` 
python pred_ld.py --file-path PRED_LD_demo.txt --r2threshold 0.8 --pop EUR --maf 0.01 --ref all_panels 
````
