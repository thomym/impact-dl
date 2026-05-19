"""
Reformat a raw GWAS summary statistics file and merge MAF from gnomAD.

INPUT
    ${gwas_root}/${discovery}/gwas_raw.tsv
        A whitespace-separated GWAS summary stats file with a header row.
        Column names may be in any common convention — see the d_header dict in
        reformat_gwas() for the recognized aliases (e.g. 'beta', 'EFFECT',
        'Beta' all map to 'BETA'; 'rs_id', 'SNPID', 'MarkerName' all map to
        'SNP'; etc.). If your file uses headers not in that dict, either rename
        them manually before running OR add the mapping to d_header.
        Required (after renaming): SNP, CHR, BP, A1, A2, P. Optional: SE, INFO,
        N, BETA, OR, MAF (filled with sensible defaults if missing).

OUTPUT
    ${gwas_root}/${discovery}/gwas.tsv
        12 columns in this order:
            SNP CHR BP A1 A2 MAF SE P N INFO BETA OR

USAGE
    python prepare_gwas.py --discovery <gwas_name>

    Paths default to gwas_root / resources / gnomad_version from paths.yaml.
    See `--help` for per-GWAS knobs and overrides.
"""

import numpy as np
import argparse
import os
import pandas as pd
from multiprocessing import Pool


def merge_maf_per_chr(df, chrom, discovery_population, maf_th, version, gnomad_dir):

    print(f"start merging maf for chr {chrom}")
    maf_file = os.path.join(
        gnomad_dir,
        f"gnomad.{version}.chr{chrom}_{discovery_population.lower()}_maf_{maf_th.replace('.', '')}",
    )
    try:
        df_maf_original = pd.read_csv(maf_file, sep='\t')
        df_maf_original['CHROM']=df_maf_original['CHROM'].apply(lambda a: int(a.replace('chr', '')) if type(a)==str else a)
        df_maf_reverse = df_maf_original.copy()
        df_maf_reverse['REF'] = df_maf_original['ALT']
        df_maf_reverse['ALT'] = df_maf_original['REF']
        df_maf_reverse['MAF'] = 1-df_maf_original['MAF']
        df_maf=pd.concat((df_maf_original, df_maf_reverse),ignore_index=True)
    except Exception as e:
        print(f"Error reading gnomAD file for chr {chrom}: {maf_file}")
        print(f"  ({type(e).__name__}) {e}")
        raise
    df_maf.columns = [f'{c}_gnomad' for c in df_maf.columns]

    df_with_maf=pd.merge(df[(df['CHR'].astype(str)==str(chrom)) | (df['CHR'].astype(str)==f'chr{chrom}') ], df_maf, how='left', left_on=['CHR', 'BP', 'A2', 'A1'], right_on=['CHROM_gnomad', 'POS_gnomad', 'REF_gnomad', 'ALT_gnomad'])

    df_with_maf = df_with_maf.apply(align_maf, axis=1)
    print(df_with_maf.head())
    res=(chrom, df_with_maf)
    print(f"end merging maf for chr {chrom}")
    return res 

def align_maf(row):
    if row['MAF_gnomad'] > 0.5:
        row_copy=row.copy()
        row['MAF_gnomad'] = 1 - row_copy['MAF_gnomad']
        row['A1'] = row_copy['A2']
        row['A2'] = row_copy['A1']
        row['BETA'] = -row_copy['BETA']
        row['OR'] = np.exp(-row_copy['BETA'])

    row['MAF']=row['MAF_gnomad']
    return row

def reformat_gwas(discovery, discovery_population, N, default_maf, maf_th, chrs, version, work_dir, gnomad_dir):
    print(f"starting reformatting {discovery}")
    df = pd.read_csv(os.path.join(work_dir, discovery, 'gwas_raw.tsv'), sep=r'\s+')
    headers = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'MAF', 'SE', 'P', 'N', 'INFO', 'BETA', 'OR']
    print(f"# of SNPs in raw GWAS: {len(df)}")

    # A1, A2
    if all([f in df.columns for f in ['A1', 'ALT', 'REF']]):
        df.loc[:, 'A2'] = df.apply(lambda a: a['REF'] if a['ALT'] == a['A1'] else a['ALT'], axis=1)
        df.drop(['REF', 'ALT'], axis=1, inplace=True)

    # Reformat A1, A2
    if all([f in df.columns for f in ['A1', 'ALT', 'REF']]):
        df.loc[:, 'A2'] = df.apply(lambda a: a['REF'] if a['ALT'] == a['A1'] else a['ALT'], axis=1)


    d_header = {}
    d_header.update({k: 'OR' for k in ['OR(A1)', 'or', 'odds_ratio']})
    d_header.update({k: 'BETA' for k in ['Beta', 'effect_size', 'beta', 'Beta' ,'EFFECT', 'Effect']})
    d_header.update({k: 'SNP' for k in ['MarkerName', 'SNPID', 'rs_id', 'snpid', 'variant_id', 'ID']})
    d_header.update({k: 'CHR' for k in ['Chr', 'Chromosome', 'chr', 'chromosome', 'CHROM', '#CHROM']})
    d_header.update({k: 'N' for k in ['sample_size', 'OBS_CT']})
    d_header.update({k: 'P' for k in ['P-val', 'Pval', 'Pvalue', 'pval', 'pvalue', 'p.value', 'P.value', 'P_BOLT_LMM', 'p_value']})
    d_header.update({k: 'BP' for k in ['POS', 'Position(hg19)', 'Position', 'pos', 'base_pair_location_grch37', 'base_pair_location']})
    d_header.update({k: 'SE' for k in ['standared', 'se', 'standard_error', 'STDERR', 'LOG(OR)_SE']})
    d_header.update({k: 'A1' for k in ['Effect_allele', 'effect_allele', 'a1', 'ALLELE1', 'ALT', 'Allele1']})
    d_header.update({k: 'A2' for k in ['Non_Effect_allele', 'noneffect_allele', 'other_allele', 'ALLELE0', 'REF', 'Allele2']})
    # d_header.update({k: 'MAF' for k in ['effect_allele_frequency']})

    df.rename(columns=d_header, inplace=True)

    # Fill missing INFO col, if necessary (default value: 1.0):
    if not 'INFO' in df.columns:
        df.loc[:, 'INFO'] = 1.0

    # Fill missing SE col, if necessary (default value: 0.005):
    if not 'SE' in df.columns:
        df.loc[:, 'SE'] = 0.005

    # Fill missing BETA col using OR col
    if not 'BETA' in df.columns and 'OR' in df.columns:
        print(df[df['OR'].apply(pd.to_numeric, errors='coerce').isnull()])
        df.loc[:, 'BETA'] = np.log(df.loc[:, 'OR'])

    # Fill missing OR col using BETA col
    if not 'OR' in df.columns and 'BETA' in df.columns:
        df.loc[:, 'OR'] = np.exp(df.loc[:, 'BETA'])

    # Fill missing N col, if necessary (default value provided by user):
    if 'Nca' in df.columns and 'Nco' in df.columns:
        df.loc[:, 'N']=df.loc[:,'Nca']+df.loc[:,'Nco']
    if not N is None and not 'N' in df.columns:        
         df.loc[:, 'N'] = N

    print(f"# of remaining SNPs: {len(df)}")
    print('start merging MAF')
    params=[]
    for chrom in chrs:
        params.append((df, chrom, discovery_population, maf_th, version, gnomad_dir))
    
    with Pool(22) as p:
        chrom_with_mafs=p.starmap(merge_maf_per_chr, params)
        print("concating chroms")
        df=pd.concat([a[1] for a in sorted(chrom_with_mafs, key=lambda a: a[0])])
        print("Done concating chroms")


    print(f'Didn\'t find matched MAF for {df.loc[:, "MAF"].isnull().sum()}/{len(df)} SNPs. Setting MAF={default_maf} as default')
    df['MAF'].fillna(default_maf, inplace=True)
    print('End merging MAF')

    df=df.reindex(headers, axis=1).replace([np.inf, -np.inf], np.nan)
    print(f'Warning: found null in {df.isnull().any(axis=1).sum()}/{len(df)} rows.')
    df.to_csv(os.path.join(work_dir, discovery, 'gwas.tsv'), index=False, sep='\t')
    print(f'Done reformatting {discovery}')


if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from paths import load_paths

    parser = argparse.ArgumentParser(description='Reformat a GWAS summary stats file and merge MAF from gnomAD.')
    parser.add_argument('-d', '--discovery', dest='discovery', required=True,
                        help='GWAS name (subdirectory under gwas_root). Comma-separated for multiple.')
    parser.add_argument('-n', '--N', dest='N', default=100000,
                        help='Sample size to fill the N column if missing.')
    parser.add_argument('-t', '--maf_th', dest='maf_th', default="0",
                        help='MAF threshold used in the gnomAD filename suffix.')
    parser.add_argument('-c', '--chrs', dest='chrs', default="all",
                        help='Chromosomes to process (e.g. "1,2,3" or "all").')
    parser.add_argument('-dp', '--discovery_population', dest='discovery_population', default='EUR',
                        help='Population code used in the gnomAD filename.')
    parser.add_argument('-dm', '--default_maf', dest='default_maf', default=0.05,
                        help='Fallback MAF when no gnomAD match is found.')
    # Optional overrides — defaults come from paths.yaml.
    parser.add_argument('-v', '--version', dest='version', default=None,
                        help='gnomAD version (default: gnomad_version from paths.yaml).')
    parser.add_argument('-w', '--work_dir', dest='work_dir', default=None,
                        help='Root containing one subdir per GWAS (default: gwas_root from paths.yaml).')
    parser.add_argument('-a', '--gnomad_dir', dest='gnomad_dir', default=None,
                        help='Directory holding per-chromosome gnomAD MAF tables '
                             '(default: gnomad_dir from paths.yaml).')
    parser.add_argument('--paths_yaml', default=None,
                        help='Path to paths.yaml (default: <repo>/paths.yaml).')
    args = parser.parse_args()

    paths = load_paths(args.paths_yaml)
    work_dir = args.work_dir or paths['gwas_root']
    gnomad_dir = args.gnomad_dir or paths['gnomad_dir']
    version = args.version or paths['gnomad_version']

    gwass = args.discovery.split(',')
    chrs = list(range(1, 23)) if args.chrs == "all" else [int(a) for a in args.chrs.split(',')]
    N = args.N
    discovery_population = args.discovery_population
    maf_th = args.maf_th
    default_maf = float(args.default_maf)
    for discovery in gwass:
        reformat_gwas(discovery, discovery_population, N, default_maf, maf_th, chrs, version, work_dir, gnomad_dir)