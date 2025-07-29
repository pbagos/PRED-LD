import dask.dataframe as dd
import numpy as np
import pandas as pd
import gc
import polars as pl

pd.set_option('display.width', None)




def TOP_LD_info_pairwise(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list=None):
    print("Loading TOP-LD files...")

    # 1) Combine rsIDs up‐front
    if imp_snp_list:
        all_rsids = list(set(rs_list) | set(imp_snp_list))
    else:
        all_rsids = rs_list

    # 2) Lazy‐scan the MAF file, project only needed cols, then filter
    maf_path = (
        f"D:/ref_parquet/ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet"
    )
    maf_lazy = (
        pl.scan_parquet(maf_path)
        .select(["Position", "rsID", "MAF", "REF", "ALT"])
        .filter(pl.col("MAF") >= maf_threshold)
        .filter(pl.col("rsID").is_in(all_rsids))
    )

    # 3) Lazy‐scan the LD file, project and filter by R²
    ld_path = (
        f"D:/ref_parquet/ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet"
    )
    ld_lazy = (
        pl.scan_parquet(ld_path)
        .select(["SNP1", "SNP2", "R2", "+/-corr", "Dprime"])
        .filter(pl.col("R2") >= R2_threshold)
    )

    # 4) Prepare two MAF views for joining
    maf1 = maf_lazy.rename({
        "Position": "SNP1", "rsID": "rsID1", "MAF": "MAF1",
        "REF": "REF1", "ALT": "ALT1"
    })
    maf2 = maf_lazy.rename({
        "Position": "SNP2", "rsID": "rsID2", "MAF": "MAF2",
        "REF": "REF2", "ALT": "ALT2"
    })

    # 5) Join LD ↔ MAF1 ↔ MAF2 all lazily
    joined = (
        ld_lazy
        .join(maf1, on="SNP1", how="inner")
        .join(maf2, on="SNP2", how="inner")
    )

    # 6) If you provided an imp_snp_list, filter rsID roles
    if imp_snp_list:
        joined = joined.filter(
            (pl.col("rsID1").is_in(rs_list)) &
            (pl.col("rsID2").is_in(imp_snp_list))
        )

    # 7) Select + rename final output columns
    final_lazy = joined.select([
        pl.col("SNP1").alias("pos1"),
        pl.col("SNP2").alias("pos2"),
        "R2", "+/-corr", "Dprime",
        "rsID1", "rsID2",
        "MAF1", "MAF2",
        "REF1", "ALT1", "REF2", "ALT2"
    ])

    # 8) Execute once in streaming mode (__very__ memory‐efficient)

    result = final_lazy.collect()

    if result.is_empty():
        print("No SNPs found.")

    # 10) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    result = result.filter(
        pl.col("rsID1").is_in(rs_list) &
        pl.col("rsID2").is_in(rs_list)
    )


    # 9) Cleanup
    del maf_lazy, maf1, maf2, ld_lazy, joined, final_lazy
    gc.collect()

    return result


def TOP_LD_process_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data

    outputData = TOP_LD_info_pairwise(list(study_df['SNP']), chromosome, population, maf_input, r2threshold, imp_snp_list)

    return outputData




def Hap_Map_LD_info_dask_pairwise(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading Hap Map files ({population}) ...")


    # 1) Validate population
    valid_pops = ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK',
                  'CHD', 'GIH', 'TSI', 'MEX', 'ASW']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    maf_file = f'D:/ref_parquet/ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.parquet'
    ld_file = f'D:/ref_parquet/ref/Hap_Map/ld_chr{chrom}_{population}.parquet'

    # 3) Read and preprocess MAF table
    maf_df = (
        pl.read_parquet(maf_file)
        # compute MAF as 'otherallele_freq'
        .with_columns(pl.col("otherallele_freq").alias("MAF"))
        # keep only the needed columns
        .select(["rs#", "refallele", "otherallele", "MAF"])
        # rename to match downstream joins
        .rename({
            "rs#": "rsID",
            "refallele": "REF",
            "otherallele": "ALT"
        })
        # filter on MAF threshold
        .filter(pl.col("MAF") >= float(maf_threshold))
    )

    # 4) Read and preprocess LD table
    #    Parquet likely comes with integer column names (0,1,2,…). We filter & rename by position.
    ld_raw = pl.read_parquet(ld_file)
    cols = ld_raw.columns  # e.g. ['0','1','2',…] or positional ints

    ld_df = (
        ld_raw
        # keep only pairs where alleles differ and R2 >= threshold
        .filter(
            (pl.col(cols[3]) != pl.col(cols[4])) &
            (pl.col(cols[6]) >= float(R2_threshold))
        )
        # rename the first 7 columns to names
        .rename({
            cols[0]: "pos1",
            cols[1]: "pos2",
            cols[2]: "pop",
            cols[3]: "rsID1",
            cols[4]: "rsID2",
            cols[5]: "Dprime",
            cols[6]: "R2"
        })
        # drop any extra columns beyond the 7 we care about
        .select(["pos1", "pos2", "pop", "rsID1", "rsID2", "Dprime", "R2"])
    )

    # 5) Join LD ↔ MAF on rsID1 and rsID2
    #    First join on rsID1 (adds REF, ALT, MAF from maf_df)
    tmp = ld_df.join(maf_df, left_on="rsID1", right_on="rsID", how="inner")
    #    Then join that result on rsID2 (suffix="_2" for the second maf columns)
    merged = tmp.join(
        maf_df,
        left_on="rsID2",
        right_on="rsID",
        how="inner",
        suffix="_2"
    )

    # 6) Select & reorder columns, and rename duplicates
    result = merged.select([
        pl.col("pos1"),
        pl.col("pos2"),
        pl.col("rsID1"),
        pl.col("rsID2"),
        pl.col("MAF").alias("MAF1"),
        pl.col("MAF_2").alias("MAF2"),
        pl.col("REF").alias("REF1"),
        pl.col("REF_2").alias("REF2"),
        pl.col("ALT").alias("ALT1"),
        pl.col("ALT_2").alias("ALT2"),
        pl.col("R2"),
        pl.col("Dprime")
    ])

    # 7) Apply the user’s SNP‐list filter
    if imp_snp_list:
        result = result.filter(
            pl.col("rsID1").is_in(rs_list) &
            pl.col("rsID2").is_in(imp_snp_list)
        )
    else:
        result = result.filter(pl.col("rsID1").is_in(rs_list))

    # # 8) Rename MAF columns to ALT_AF and convert to Pandas
    # final_pl = result.rename({
    # "MAF1": "ALT_AF1",
    # "MAF2": "ALT_AF2"
    # })

    # Convert to Pandas and reset index
    # 10) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    result = result.filter(
        pl.col("rsID1").is_in(rs_list) &
        pl.col("rsID2").is_in(rs_list)
    )

    # 9) Clean up
    del ld_raw, ld_df, maf_df, tmp, merged
    gc.collect()

    # if result.empty:
    #     print("No SNPs found")




    return result


def Hap_Map_process_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = Hap_Map_LD_info_dask_pairwise(list(study_df['SNP']), chromosome, population, maf_input, r2threshold,
                                      imp_snp_list)



    return outputData


def pheno_Scanner_LD_info_dask_pairwise(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    # 1) Enforce minimum R2 threshold
    if R2_threshold < 0.8:
        print("Pheno Scanner includes data with R2 ≥ 0.8. Setting R2_threshold = 0.8")
        R2_threshold = 0.8

    print("Building lazy plan for Pheno Scanner files...")

    # 2) Paths & pop mapping
    maf_file = "D:/ref_parquet/ref/Pheno_Scanner/1000G.parquet"
    ld_file = f"D:/ref_parquet/ref/Pheno_Scanner/1000G_{population}/1000G_{population}_chr{chrom}.parquet"
    pop_map = {"EUR": "eur", "EAS": "eas", "AFR": "afr", "AMR": "amr", "SAS": "sas"}
    maf_pop = pop_map.get(population)
    if maf_pop is None:
        raise ValueError(f"Unsupported population: {population}")

    # 3) First MAF scan (for ref side → MAF1)
    maf1_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter(
            (pl.col("chr") == chrom) &
            (pl.col(maf_pop) != "-")
        )
        .with_columns(pl.col(maf_pop).cast(pl.Float64))
        .filter(pl.col(maf_pop) >= maf_threshold)
        .rename({
            "hg19_coordinates": "ref_hg19_coordinates",
            "rsid": "ref_rsid",
            maf_pop: "MAF1",
            "a1": "ALT1",
            "a2": "REF1",
        })
        .select(["ref_hg19_coordinates", "ref_rsid", "MAF1", "ALT1", "REF1"])
    )

    # 4) LD scan
    ld_lazy = (
        pl.scan_parquet(ld_file)
        .select(["ref_hg19_coordinates", "ref_rsid", "rsid", "r2", "r", "dprime"])
        .filter(
            (pl.col("ref_rsid") != pl.col("rsid")) &
            (pl.col("r2") >= R2_threshold)
        )
    )

    # 5) Second MAF scan (for query side → MAF2), keep `rsid` for join!
    maf2_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter(
            (pl.col("chr") == chrom) &
            (pl.col(maf_pop) != "-")
        )
        .with_columns(pl.col(maf_pop).cast(pl.Float64))
        .filter(pl.col(maf_pop) >= maf_threshold)
        .rename({
            "hg19_coordinates": "pos2(hg19)",
            maf_pop: "MAF2",
            "a1": "ALT2",
            "a2": "REF2",
        })
        # note: we keep `rsid` here so we can join on it
        .select(["pos2(hg19)", "rsid", "MAF2", "ALT2", "REF2"])
    )

    # 6) Join the three pieces
    joined = (
        ld_lazy
        .join(maf1_lazy, on="ref_rsid", how="inner")
        .join(maf2_lazy, on="rsid", how="inner")
    )

    # 7) Rename LD columns & the `rsid` from maf2 → final names
    renamed = joined.rename({
        "ref_rsid": "rsID1",
        "rsid": "rsID2",
        "ref_hg19_coordinates": "pos1(hg19)",
        "r2": "R2",
    })

    # 8) Filter by user‐supplied SNP lists
    filtered = renamed.filter(
        pl.col("rsID2").is_in(rs_list) &
        (pl.col("rsID1").is_in(imp_snp_list) if imp_snp_list else pl.lit(True))
    )

    # 9) Extract numeric coords
    with_pos = filtered.with_columns([
        pl.col("pos1(hg19)").str.extract(r".*:(\d+)$", 1).cast(pl.Int64),
        pl.col("pos2(hg19)").str.extract(r".*:(\d+)$", 1).cast(pl.Int64),
    ])

    # 10) Final projection/order
    final_lazy = with_pos.select([
        "rsID1", "pos1(hg19)", "rsID2", "dprime",
        "pos2(hg19)", "R2", "r",
        "MAF1", "MAF2", "ALT1", "REF1", "ALT2", "REF2",
    ])

    # 11) Execute in streaming, collect to Polars → pandas
    final_pl = final_lazy.collect()


    final_pl = final_pl.filter(
        pl.col("rsID1").is_in(rs_list) &
        pl.col("rsID2").is_in(rs_list)
    )

    # 12) Cleanup
    gc.collect()
    if final_pl.empty:
        print("No SNPs found")

    return final_pl


def pheno_Scanner_process_pairwise(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = pheno_Scanner_LD_info_dask_pairwise(list(study_df['SNP']), chromosome, population, maf_input, r2threshold,
                                            imp_snp_list)

    return outputData



def hg38_1kg_LD_info_pairwise(
    rs_list: list[str],
    chrom: int,
    population: str,
    maf_threshold: float | None,
    R2_threshold: float | None,
    imp_snp_list: list[str] | None
) -> pl.DataFrame:
    """
    Fetch pairwise LD (R and R2) for all variant pairs on `chrom` in `population`
    where both SNPs are in `rs_list` (or in `imp_snp_list` if provided),
    apply MAF and R2 thresholds, and return a Polars DataFrame with allele info.
    """
    print(f"Loading 1000 Genomes Project (hg38) files ({population}) ...")

    # 1) Validate population
    valid_pops = ['AMR', 'EUR', 'AFR', 'SAS', 'EAS']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    maf_file = f'D:/ref_parquet/ref/1000G_hg38/1000G_{population}_0_01.parquet'
    ld_file  = f'D:/ref_parquet/ref/1000G_hg38/{population}/chr{chrom}_merged.parquet'

    # 3) Load & filter MAF data
    print("Loading and filtering MAF file...")
    tmp = pl.read_parquet(maf_file)
    maf_df = tmp.rename({
        old: new for old, new in zip(tmp.columns[:6],
                                     ['CHR','SNP','MAF','POS','REF','ALT'])
    })
    if maf_threshold is not None:
        maf_df = maf_df.filter(pl.col("MAF") >= maf_threshold)

    # decide which SNP‐set to use
    snp_set = set(imp_snp_list) if imp_snp_list else set(rs_list)
    maf_df = maf_df.filter(pl.col("SNP").is_in(snp_set))

    # 4) Load & filter LD data
    print("Loading and filtering LD data...")
    ld_df = (
        pl.read_parquet(ld_file, columns=["CHR_A","SNP_A","CHR_B","SNP_B","R"])
          .with_columns((pl.col("R") ** 2).alias("R2"))
          # keep only pairs where both sides are in your SNP set
          .filter(
              pl.col("SNP_A").is_in(snp_set) &
              pl.col("SNP_B").is_in(snp_set)
          )
    )
    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)

    # 5) Annotate allele A
    merged = (
        ld_df.join(
            maf_df,
            left_on=["CHR_A","SNP_A"],
            right_on=["CHR","SNP"],
            how="left"
        )
        .rename({
            "CHR_A":"CHR",
            "MAF":"MAF_A",
            "POS":"POS_A",
            "REF":"REF_A",
            "ALT":"ALT_A"
        })
    )

    # 6) Annotate allele B
    merged = (
        merged.join(
            maf_df,
            left_on=["CHR_B","SNP_B"],
            right_on=["CHR","SNP"],
            how="left"
        )
        .rename({
            "MAF":"MAF_B",
            "POS":"POS_B",
            "REF":"REF_B",
            "ALT":"ALT_B"
        })
        .drop(["CHR_B"])
    )

    # 7) Drop any pairs where one allele lacked MAF info
    merged = merged.drop_nulls(subset=["MAF_A","MAF_B"])

    # 8) (Re‐)apply MAF threshold on both alleles, if requested
    if maf_threshold is not None:
        merged = merged.filter(
            (pl.col("MAF_A") >= maf_threshold) &
            (pl.col("MAF_B") >= maf_threshold)
        )

    # 9) Reorder columns for readability
    preferred = [
        "SNP_A","CHR","POS_A","SNP_B","POS_B",
        "REF_A","ALT_A","MAF_A","REF_B","ALT_B","MAF_B",
        "R","R2"
    ]
    cols = merged.columns
    final_cols = [c for c in preferred if c in cols] + [c for c in cols if c not in preferred]

    # 10) FINAL FILTER: ensure both SNP_A and SNP_B are in the original snp_set
    merged = merged.filter(
        pl.col("SNP_A").is_in(snp_set) &
        pl.col("SNP_B").is_in(snp_set)
    )
   # print (merged)
    # return with your preferred column order
    return merged.select(final_cols)


def hg38_1kg_process_pairwise(
    study_df: pl.DataFrame,
    r2threshold: float | None,
    population: str,
    maf_input: float | None,
    chromosome: int,
    imp_snp_list: list[str] | None = None
) -> pl.DataFrame:
    """
    Wrapper to fetch pairwise LD for SNPs in `study_df['SNP']` on `chromosome`.
    """
    rs_list = list(study_df['SNP'])
    return hg38_1kg_LD_info_pairwise(
        rs_list=rs_list,
        chrom=chromosome,
        population=population,
        maf_threshold=maf_input,
        R2_threshold=r2threshold,
        imp_snp_list=imp_snp_list
    )





def hg38_1kg_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading 1000 Genomes Project (hg38) files ({population}) ...")

    # 1) Validate population
    valid_pops = ['AMR', 'EUR', 'AFR', 'SAS', 'EAS']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths

    maf_file = f'D:/ref_parquet/ref/1000G_hg38/1000G_{population}_0_01.parquet'
    ld_file = f'D:/ref_parquet/ref/1000G_hg38/{population}/chr{chrom}_merged.parquet'

     # 3) Load and filter MAF file early
    print("Loading and filtering MAF file...")



    maf_df = pl.read_parquet(maf_file ).rename({
        maf_col: new_col for maf_col, new_col in zip(
            pl.read_parquet(maf_file).columns[:6],
            ['CHR', 'SNP', 'MAF', 'POS', 'REF', 'ALT']
        )
    })

    if maf_threshold is not None:
        maf_df = maf_df.filter(pl.col("MAF") >= maf_threshold)

    if imp_snp_list:
        maf_df = maf_df.filter(pl.col("SNP").is_in(imp_snp_list))
    elif rs_list:
        maf_df = maf_df.filter(pl.col("SNP").is_in(rs_list))

    maf_snps = maf_df.select("SNP").to_series().to_list()

    # 4) Load LD data
    print("Loading and filtering LD data...")
    ld_df = pl.read_parquet(ld_file, columns=["CHR_A", "SNP_A", "CHR_B", "SNP_B", "R"])

    if imp_snp_list or rs_list:
        filter_snps = set(imp_snp_list or rs_list)
        ld_df = ld_df.filter(
            pl.col("SNP_A").is_in(filter_snps) | pl.col("SNP_B").is_in(filter_snps)
        )

    ld_df = ld_df.with_columns((pl.col("R") ** 2).alias("R2"))

    if R2_threshold is not None:
        ld_df = ld_df.filter(pl.col("R2") >= R2_threshold)

    # 5) Join with allele A
    merged = ld_df.join(
        maf_df,
        left_on=["CHR_A", "SNP_A"],
        right_on=["CHR", "SNP"],
        how="left"
    ).rename({
        "CHR_A": "CHR",
        "MAF": "MAF_A",
        "POS": "POS_A",
        "REF": "REF_A",
        "ALT": "ALT_A"
    }) # Remove duplicated join columns

    # 6) Join with allele B
    merged = merged.join(
        maf_df,
        left_on=["CHR_B", "SNP_B"],
        right_on=["CHR", "SNP"],
        how="left"
    ).rename({
        "MAF": "MAF_B",
        "POS": "POS_B",
        "REF": "REF_B",
        "ALT": "ALT_B"
    }).drop(["CHR_B"])  # Drop redundant CHR and SNP columns

    # 7) Final MAF filtering
    if maf_threshold is not None:
        merged = merged.filter(
            (pl.col("MAF_A") >= maf_threshold) & (pl.col("MAF_B") >= maf_threshold)
        )

    # 8) Reorder columns
    preferred_order = [
        "SNP_A", "CHR", "POS_A", "SNP_B", "POS_B",
        "REF_A", "ALT_A", "MAF_A", "REF_B", "ALT_B", "MAF_B",
        "R", "R2"
    ]
    # Ensure all columns exist before reordering
    existing_columns = merged.columns
    final_order = [col for col in preferred_order if col in existing_columns] + [
        col for col in existing_columns if col not in preferred_order
    ]
    merged = merged.select(final_order)

    return merged


def hg38_1kg_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data

    outputData = hg38_1kg_LD_info(list(study_df['SNP']), chromosome, population, maf_input, r2threshold, imp_snp_list)


    return outputData

def Hap_Map_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    print(f"Loading Hap Map files ({population}) ...")


    # 1) Validate population
    valid_pops = ['YRI', 'CHB', 'JPT', 'CEU', 'MKK', 'LWK',
                  'CHD', 'GIH', 'TSI', 'MEX', 'ASW']
    if population not in valid_pops:
        raise ValueError(f"Population '{population}' not available. Choose from: {valid_pops}")

    # 2) File paths
    maf_file = f'D:/ref_parquet/ref/Hap_Map/allele_freqs_chr{chrom}_{population}_phase3.2_nr.b36_fwd.parquet'
    ld_file = f'D:/ref_parquet/ref/Hap_Map/ld_chr{chrom}_{population}.parquet'

    # 3) Read and preprocess MAF table
    maf_df = (
        pl.read_parquet(maf_file)
        # compute MAF as 'otherallele_freq'
        .with_columns(pl.col("otherallele_freq").alias("MAF"))
        # keep only the needed columns
        .select(["rs#", "refallele", "otherallele", "MAF"])
        # rename to match downstream joins
        .rename({
            "rs#": "rsID",
            "refallele": "REF",
            "otherallele": "ALT"
        })
        # filter on MAF threshold
        .filter(pl.col("MAF") >= float(maf_threshold))
    )

    # 4) Read and preprocess LD table
    #    Parquet likely comes with integer column names (0,1,2,…). We filter & rename by position.
    ld_raw = pl.read_parquet(ld_file)
    cols = ld_raw.columns  # e.g. ['0','1','2',…] or positional ints

    ld_df = (
        ld_raw
        # keep only pairs where alleles differ and R2 >= threshold
        .filter(
            (pl.col(cols[3]) != pl.col(cols[4])) &
            (pl.col(cols[6]) >= float(R2_threshold))
        )
        # rename the first 7 columns to names
        .rename({
            cols[0]: "pos1",
            cols[1]: "pos2",
            cols[2]: "pop",
            cols[3]: "rsID1",
            cols[4]: "rsID2",
            cols[5]: "Dprime",
            cols[6]: "R2"
        })
        # drop any extra columns beyond the 7 we care about
        .select(["pos1", "pos2", "pop", "rsID1", "rsID2", "Dprime", "R2"])
    )

    # 5) Join LD ↔ MAF on rsID1 and rsID2
    #    First join on rsID1 (adds REF, ALT, MAF from maf_df)
    tmp = ld_df.join(maf_df, left_on="rsID1", right_on="rsID", how="inner")
    #    Then join that result on rsID2 (suffix="_2" for the second maf columns)
    merged = tmp.join(
        maf_df,
        left_on="rsID2",
        right_on="rsID",
        how="inner",
        suffix="_2"
    )

    # 6) Select & reorder columns, and rename duplicates
    result = merged.select([
        pl.col("pos1"),
        pl.col("pos2"),
        pl.col("rsID1"),
        pl.col("rsID2"),
        pl.col("MAF").alias("MAF1"),
        pl.col("MAF_2").alias("MAF2"),
        pl.col("REF").alias("REF1"),
        pl.col("REF_2").alias("REF2"),
        pl.col("ALT").alias("ALT1"),
        pl.col("ALT_2").alias("ALT2"),
        pl.col("R2"),
        pl.col("Dprime")
    ])

    # 7) Apply the user’s SNP‐list filter
    if imp_snp_list:
        result = result.filter(
            pl.col("rsID1").is_in(rs_list) &
            pl.col("rsID2").is_in(imp_snp_list)
        )
    else:
        result = result.filter(pl.col("rsID1").is_in(rs_list))

    # # 8) Rename MAF columns to ALT_AF and convert to Pandas
    # final_pl = result.rename({
    # "MAF1": "ALT_AF1",
    # "MAF2": "ALT_AF2"
    # })

    # # Convert to Pandas and reset index
    # final_pd = result.to_pandas().reset_index(drop=True)

    # 9) Clean up
    del ld_raw, ld_df, maf_df, tmp, merged
    gc.collect()



    return result


def Hap_Map_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = Hap_Map_LD_info_dask(list(study_df['SNP']), chromosome, population, maf_input, r2threshold,
                                      imp_snp_list)



    return outputData


def pheno_Scanner_LD_info_dask(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list):
    # 1) Enforce minimum R2 threshold
    if R2_threshold < 0.8:
        print("Pheno Scanner includes data with R2 ≥ 0.8. Setting R2_threshold = 0.8")
        R2_threshold = 0.8

    print("Building lazy plan for Pheno Scanner files...")

    # 2) Paths & pop mapping
    maf_file = "D:/ref_parquet/ref/Pheno_Scanner/1000G.parquet"
    ld_file = f"D:/ref_parquet/ref/Pheno_Scanner/1000G_{population}/1000G_{population}_chr{chrom}.parquet"
    pop_map = {"EUR": "eur", "EAS": "eas", "AFR": "afr", "AMR": "amr", "SAS": "sas"}
    maf_pop = pop_map.get(population)
    if maf_pop is None:
        raise ValueError(f"Unsupported population: {population}")

    # 3) First MAF scan (for ref side → MAF1)
    maf1_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter(
            (pl.col("chr") == chrom) &
            (pl.col(maf_pop) != "-")
        )
        .with_columns(pl.col(maf_pop).cast(pl.Float64))
        .filter(pl.col(maf_pop) >= maf_threshold)
        .rename({
            "hg19_coordinates": "ref_hg19_coordinates",
            "rsid": "ref_rsid",
            maf_pop: "MAF1",
            "a1": "ALT1",
            "a2": "REF1",
        })
        .select(["ref_hg19_coordinates", "ref_rsid", "MAF1", "ALT1", "REF1"])
    )

    # 4) LD scan
    ld_lazy = (
        pl.scan_parquet(ld_file)
        .select(["ref_hg19_coordinates", "ref_rsid", "rsid", "r2", "r", "dprime"])
        .filter(
            (pl.col("ref_rsid") != pl.col("rsid")) &
            (pl.col("r2") >= R2_threshold)
        )
    )

    # 5) Second MAF scan (for query side → MAF2), keep `rsid` for join!
    maf2_lazy = (
        pl.scan_parquet(maf_file)
        .select(["hg19_coordinates", "chr", "rsid", maf_pop, "a1", "a2"])
        .filter(
            (pl.col("chr") == chrom) &
            (pl.col(maf_pop) != "-")
        )
        .with_columns(pl.col(maf_pop).cast(pl.Float64))
        .filter(pl.col(maf_pop) >= maf_threshold)
        .rename({
            "hg19_coordinates": "pos2(hg19)",
            maf_pop: "MAF2",
            "a1": "ALT2",
            "a2": "REF2",
        })
        # note: we keep `rsid` here so we can join on it
        .select(["pos2(hg19)", "rsid", "MAF2", "ALT2", "REF2"])
    )

    # 6) Join the three pieces
    joined = (
        ld_lazy
        .join(maf1_lazy, on="ref_rsid", how="inner")
        .join(maf2_lazy, on="rsid", how="inner")
    )

    # 7) Rename LD columns & the `rsid` from maf2 → final names
    renamed = joined.rename({
        "ref_rsid": "rsID1",
        "rsid": "rsID2",
        "ref_hg19_coordinates": "pos1(hg19)",
        "r2": "R2",
    })

    # 8) Filter by user‐supplied SNP lists
    filtered = renamed.filter(
        pl.col("rsID2").is_in(rs_list) &
        (pl.col("rsID1").is_in(imp_snp_list) if imp_snp_list else pl.lit(True))
    )

    # 9) Extract numeric coords
    with_pos = filtered.with_columns([
        pl.col("pos1(hg19)").str.extract(r".*:(\d+)$", 1).cast(pl.Int64),
        pl.col("pos2(hg19)").str.extract(r".*:(\d+)$", 1).cast(pl.Int64),
    ])

    # 10) Final projection/order
    final_lazy = with_pos.select([
        "rsID1", "pos1(hg19)", "rsID2", "dprime",
        "pos2(hg19)", "R2", "r",
        "MAF1", "MAF2", "ALT1", "REF1", "ALT2", "REF2",
    ])

    # 11) Execute in streaming, collect to Polars → pandas
    final_pl = final_lazy.collect()
    final_pd = final_pl.to_pandas()

    # 12) Cleanup
    gc.collect()
    if final_pd.empty:
        print("No SNPs found")

    return final_pd


def pheno_Scanner_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data
    outputData = pheno_Scanner_LD_info_dask(list(study_df['SNP']), chromosome, population, maf_input, r2threshold,
                                            imp_snp_list)

    return outputData


def TOP_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list=None):
    print("Loading TOP-LD files...")

    # 1) Combine rsIDs up‐front
    if imp_snp_list:
        all_rsids = list(set(rs_list) | set(imp_snp_list))
    else:
        all_rsids = rs_list

    # 2) Lazy‐scan the MAF file, project only needed cols, then filter
    maf_path = (
        f"D:/ref_parquet/ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet"
    )
    maf_lazy = (
        pl.scan_parquet(maf_path)
        .select(["Position", "rsID", "MAF", "REF", "ALT"])
        .filter(pl.col("MAF") >= maf_threshold)
        .filter(pl.col("rsID").is_in(all_rsids))
    )

    # 3) Lazy‐scan the LD file, project and filter by R²
    ld_path = (
        f"D:/ref_parquet/ref/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet"
    )
    ld_lazy = (
        pl.scan_parquet(ld_path)
        .select(["SNP1", "SNP2", "R2", "+/-corr", "Dprime"])
        .filter(pl.col("R2") >= R2_threshold)
    )

    # 4) Prepare two MAF views for joining
    maf1 = maf_lazy.rename({
        "Position": "SNP1", "rsID": "rsID1", "MAF": "MAF1",
        "REF": "REF1", "ALT": "ALT1"
    })
    maf2 = maf_lazy.rename({
        "Position": "SNP2", "rsID": "rsID2", "MAF": "MAF2",
        "REF": "REF2", "ALT": "ALT2"
    })

    # 5) Join LD ↔ MAF1 ↔ MAF2 all lazily
    joined = (
        ld_lazy
        .join(maf1, on="SNP1", how="inner")
        .join(maf2, on="SNP2", how="inner")
    )

    # 6) If you provided an imp_snp_list, filter rsID roles
    if imp_snp_list:
        joined = joined.filter(
            (pl.col("rsID1").is_in(rs_list)) &
            (pl.col("rsID2").is_in(imp_snp_list))
        )

    # 7) Select + rename final output columns
    final_lazy = joined.select([
        pl.col("SNP1").alias("pos1"),
        pl.col("SNP2").alias("pos2"),
        "R2", "+/-corr", "Dprime",
        "rsID1", "rsID2",
        "MAF1", "MAF2",
        "REF1", "ALT1", "REF2", "ALT2"
    ])

    # 8) Execute once in streaming mode (__very__ memory‐efficient)

    result = final_lazy.collect()

    if result.is_empty():
        print("No SNPs found.")

    # 9) Cleanup
    del maf_lazy, maf1, maf2, ld_lazy, joined, final_lazy
    gc.collect()

    return result.to_pandas()


def TOP_LD_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list):
    # Fetch LD info data

    outputData = TOP_LD_info(list(study_df['SNP']), chromosome, population, maf_input, r2threshold, imp_snp_list)

    return outputData


def process_data(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list):
    final_results_list = []

    study_df = pd.read_csv(file_path, sep="\t")
    print(study_df)
    chroms = list(set(study_df['CHR']))

    # Take only 22 chromosomes
    chroms = [chrom for chrom in chroms if not str(chrom).isdigit() or int(chrom) <= 22]
    ref_panel = ref_file




    # Depending on the reference panel...
    if ref_panel == '1000G_hg38':
        for chrom in chroms:
            final_data = hg38_1kg_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            # Save individual chromosome LD info if needed
            # final_data.write_csv(f"LD_info_chr{chrom}.txt", separator="\t")
            # print(f"Check 'LD_info_chr{chrom}.txt' for LD information")

            final_results_list.append(final_data)


        # Concatenate Polars DataFrames
        final_df = pl.concat(final_results_list)

        # Write final concatenated DataFrame to file
        final_df.write_csv("LD_info_chr_all.txt", separator="\t")
        print("Check 'LD_info_chr_all.txt' for results")

    if ref_panel == 'TOP_LD':
        for chrom in chroms:
            final_data = TOP_LD_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            # Save individual chromosome LD info if needed


            final_results_list.append(final_data)


        # Concatenate Polars DataFrames
        final_df = pd.concat(final_results_list)

        # Write final concatenated DataFrame to file
        final_df.to_csv("LD_info_chr_all.txt", sep="\t")
        print("Check 'LD_info_chr_all.txt' for results")

    if ref_panel == 'Pheno_Scanner':
        for chrom in chroms:
            final_data = pheno_Scanner_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            # Save individual chromosome LD info if needed
            final_data.to_csv(f"LD_info_chr{chrom}.txt", sep="\t")
            print(f"Check 'LD_info_chr{chrom}.txt' for LD information")

            final_results_list.append(final_data)


        final_df = pd.concat(final_results_list)

        # Write final concatenated DataFrame to file
        final_df.to_csv("LD_info_chr_all.txt", sep="\t")
        print("Check 'LD_info_chr_all.txt' for results")

    if ref_panel == 'Hap_Map':
        for chrom in chroms:
            final_data = Hap_Map_process(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            # # Save individual chromosome LD info if needed
            #
            # print(f"Check 'LD_info_chr{chrom}.txt' for LD information")

            final_results_list.append(final_data)

    
        # Concatenate Polars DataFrames
        final_df = pl.concat(final_results_list)

        # Write final concatenated DataFrame to file
        final_df.write_csv("LD_info_chr_all.txt", separator="\t")
        print("Check 'LD_info_chr_all.txt' for results")



def process_data_pairwise(file_path, r2threshold, population, maf_input, ref_file, imp_snp_list):
    final_results_list = []

    study_df = pd.read_csv(file_path, sep="\t")
    print(study_df)

    #chroms = list(set(study_df['CHR']))

    # Take only 22 chromosomes
    #  chroms = [chrom for chrom in chroms if not str(chrom).isdigit() or int(chrom) <= 22]
    ref_panel = ref_file




    # Depending on the reference panel...
    if ref_panel == '1000G_hg38':
        for chrom in range(1,23):
            final_data = hg38_1kg_process_pairwise(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            # # Save individual chromosome LD info if needed
            # final_data.write_csv(f"LD_info_chr{chrom}.txt", separator="\t")
            # print(f"Check 'LD_info_chr{chrom}.txt' for LD information")

            final_results_list.append(final_data)


        # Write final concatenated DataFrame to file
        final_df = pl.concat(final_results_list)
        final_df.write_csv("LD_info_chr_all_pairwise.txt", separator="\t")
        print("Check 'LD_info_chr_all_pairwise.txt' for results")


    if ref_panel == 'TOP_LD':
        for chrom in range(1,23):
            final_data = TOP_LD_process_pairwise(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            # # Save individual chromosome LD info if needed
            # final_data.write_csv(f"LD_info_chr{chrom}.txt", separator="\t")
            # print(f"Check 'LD_info_chr{chrom}.txt' for LD information")

            final_results_list.append(final_data)


        # Write final concatenated DataFrame to file
        final_df = pl.concat(final_results_list)
        final_df.write_csv("LD_info_chr_all_pairwise.txt", separator="\t")
        print("Check 'LD_info_chr_all_pairwise.txt' for results")

    if ref_panel == 'Hap_Map':
        for chrom in range(1,23):
            final_data = Hap_Map_process_pairwise(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            # # Save individual chromosome LD info if needed
            # final_data.write_csv(f"LD_info_chr{chrom}.txt", separator="\t")
            # print(f"Check 'LD_info_chr{chrom}.txt' for LD information")

            final_results_list.append(final_data)


        # Write final concatenated DataFrame to file
        final_df = pl.concat(final_results_list)
        final_df.write_csv("LD_info_chr_all_pairwise.txt", separator="\t")
        print("Check 'LD_info_chr_all_pairwise.txt' for results")

    if ref_panel == 'Pheno_Scanner':
        for chrom in range(1, 23):
            final_data = pheno_Scanner_process_pairwise(study_df, r2threshold, population, maf_input, chrom, imp_snp_list)

            # # Save individual chromosome LD info if needed
            # final_data.write_csv(f"LD_info_chr{chrom}.txt", separator="\t")
            # print(f"Check 'LD_info_chr{chrom}.txt' for LD information")

            final_results_list.append(final_data)

        # Write final concatenated DataFrame to file
        final_df = pl.concat(final_results_list)
        final_df.write_csv("LD_info_chr_all_pairwise.txt", separator="\t")
        print("Check 'LD_info_chr_all_pairwise.txt' for results")



