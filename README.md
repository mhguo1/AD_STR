# AD_STR

Repository for code used in Guo et al. analysis of short tandem repeats in Alzheimer's disease.

STR genotyping is performed using gangSTR v2.4.0 (available at https://github.com/gymreklab/GangSTR) and ExpansionHunter v5 (https://github.com/Illumina/ExpansionHunter). Both software programs were run on Amazon Web Services. The approximate run time using 4 Gb of RAM on AWS is 8 hours per sample for each software program.

Custom code for processing STR genotype data was performed using R v3.6.3. The software depends on the following publicly available packages in R. Approximate installation time is < 5 minutes on a MacBook Pro.
- regioneR v1.30.0
- clusterProfiler v4.2.2
- GenomicRanges v1.38.0
- dbscan v1.1-11
- RNOmni v1.0.1
- stringr v1.5.0
- reshape2 v1.4.4
- data.table v1.13.0
- RNOmni v1.0.1


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

C) Fisher's exact test for single STR burden analysis (fisher_single_str_test.R)

D) STR expansion analysis
  1) Identification of STR expansions using DBSCAN algorithm (run_dbscan.R)
  2) Parse DBSCAN results (parse_dbscan.R)

E) Genomic annotations and enrichments
 1) Perform genomic annotations
 2) Run chromHMM enrichments (run_chromhmm_enrichment.R)
 3) Run transposable element enrichments (run_te_enrichment.R)
 4) Run Hippocampus histone ChIP enrichments from ENCODE (run_encode_histone_chip_enrichment.R)


We also provide a simulated dataset for trialing the DBSCAN algorithm. This is provided along with a README in the test_dbscan folder.
