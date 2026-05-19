#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/parse_args.sh" "$@"

# Parse input
# gwas_name   -- (required) name of this GWAS; data lives in ${work_dir}/${gwas_name}/
# work_dir    -- (optional) root containing GWAS subdirectories. Defaults to gwas_root from paths.yaml.
# paths_yaml  -- (optional) override the paths.yaml location.
if [[ -z "${gwas_name}" ]]; then echo "ERROR: --gwas_name is required" >&2; exit 1; fi

if [[ -z "${work_dir}" ]]; then
    paths_args=""
    if [[ -n "${paths_yaml}" ]]; then paths_args="--paths_yaml ${paths_yaml}"; fi
    work_dir=$(python "$REPO_ROOT/paths.py" $paths_args --get gwas_root)
fi

GWAS_path="${work_dir}/${gwas_name}/"

# Run pipeline
echo filter low-quality SNPs
cat ${GWAS_path}gwas.tsv |\
awk 'NR==1 || ($6 > 0.01) && ($8 < 1.0) && ($10 > 0.8) {print}' |\
gzip  > ${GWAS_path}gwas.quality.gz

echo filter ambiguous SNPS
gunzip -c ${GWAS_path}gwas.quality.gz |\
awk '!( ($4=="A" && $5=="T") || \
		($4=="T" && $5=="A") || \
		($4=="G" && $5=="C") || \
		($4=="C" && $5=="G")) {print}' |\
	gzip > ${GWAS_path}gwas.noambig.gz

echo get duplicated SNPs
echo "dups" > ${GWAS_path}duplicated.snp
gunzip -c ${GWAS_path}gwas.noambig.gz |\
awk '{ print $1}' |\
sort |\
uniq -d >> ${GWAS_path}duplicated.snp

echo remove duplicated SNPs
awk '{if(NR==FNR) {c[$1]++; next;} if (c[$1]==0){print $0}}' <(cat ${GWAS_path}duplicated.snp) <(gunzip -c ${GWAS_path}gwas.noambig.gz) |\
gzip - > ${GWAS_path}gwas.QC.gz

# awk  grep -vf ${GWAS_path}duplicated.snp |\

gunzip -c ${GWAS_path}gwas.QC.gz > ${GWAS_path}gwas.QC.Transformed

echo extract SNP p-value
awk '{print $1,$8}' ${GWAS_path}gwas.QC.Transformed > ${GWAS_path}SNP.pvalue