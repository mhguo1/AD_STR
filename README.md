# AD_STR

Repository for code used in Guo et al. analysis of short tandem repeats in Alzheimer's disease.

Here is a summary of the code:
A) STR_genotyping
 1) ExpansionHunter genotyping (run_EH.sh)
 2) gangSTR genotying (run_gangSTR.sh)
 3) ExpnasionHunter json file (eh.v5_w_gangstr.v13.polymorphic.json.gz)
 4) gangSTR STR bed file (gangstr.v13.polymorphic_w_eh.v5_offtarget.exp10_p05.top5_no_segdup.fast_2s.bed.gz)
 5) Code to extract STR genotypes from gangSTR vcf files (parse_gangstr_genotypes.R)
 6) Code to extract coverage from gangSTR vcf files (parse_gangstr_coverage.R)

B) Single STR association analyses
  1) Run single STR associations

C) Fisher's exact test for single STR burden analysis

D) STR expansion analysis
  1) Identification of STR expansions using DBSCAN algorithm
  2) 
