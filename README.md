# AD_STR

Repository for code used in Guo et al. analysis of short tandem repeats in Alzheimer's disease.

Here is a summary of the code:

A) STR genotyping
 1) ExpansionHunter genotyping (run_EH.sh)
 2) gangSTR genotying (run_gangSTR.sh)
 3) ExpnasionHunter json file (eh.v5_w_gangstr.v13.polymorphic.json.gz)
 4) gangSTR STR bed file (gangstr.v13.polymorphic_w_eh.v5_offtarget.exp10_p05.top5_no_segdup.fast_2s.bed.gz)
 5) Code to extract coverage from ExpansionHunter vcf files (parse_eh_coverage.R)
 6) Code to extract STR genotypes from ExpansionHunter vcf files (parse_eh_genotypes.R)
 7) Code to extract STR genotypes from gangSTR vcf files (parse_gangstr_genotypes.R)
 8) Code to extract coverage from gangSTR vcf files (parse_gangstr_coverage.R)

B) Single STR association analyses (single_str_association_analysis.R)

C) Fisher's exact test for single STR burden analysis

D) STR expansion analysis
  1) Identification of STR expansions using DBSCAN algorithm (run_dbscan.R)
  2) Parse DBSCAN results (parse_dbscan.R)

E) Genomic annotations and enrichments
 1) Perform genomic annotations
 2) Run chromHMM enrichments (run_chromhmm_enrichment.R)
 3) Run transposable element enrichments (run_te_enrichment.R)
 4) Run Hippocampus histone ChIP enrichments from ENCODE (run_encode_histone_chip_enrichment.R)
