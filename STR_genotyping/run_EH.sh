#!/usr/bin/sh

#This code runs ExpansionHunter v5

#For ExpansionHunter:
ref="/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
bam=$1 #path to bam/cram
sample=$2 #sample prefix/ID
json="eh.v5_w_gangstr.v13.polymorphic.json"

/project/jcreminslab/guomic_projects/software/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter \
--reads ${bam} \
--reference ${ref} \
--variant-catalog ${json} \
--output-prefix ${sample}.eh_test \
--analysis-mode streaming \
--log-level warn
