#!/usr/bin/sh

ref="/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa" #Reference genome

str_ref="/project/jcreminslab/guomic_projects/str_test/str_ref_merge/gangstr.v13.polymorphic_w_eh.v5_offtarget.exp10_p05.top5_no_segdup.fast_2s.bed"
sample_id=$1 #sample prefixl/ID
sex=$2 #Provide sex as "F" or "M"
bam=$3 #Path to cram or bam file
samtools index ${bam}

/software/bin/GangSTR \
--bam ${bam} \
--ref ${ref} \
--regions ${str_ref} \
--out ${sample_id}.gangstr.output \
--bam-samps ${sample_id} \
--max-proc-read 100000
