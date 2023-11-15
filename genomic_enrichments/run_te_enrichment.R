#Run enrichment analysis for transposable elements
#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(SummarizedExperiment)
library(regioneR)
library(data.table)
library(stringr)


#Remove samples with high numbers of expansions
dat.dbscan.counts<-read.delim("/dbscan/eh_v5.dbscan_eps2.results_w_alleles_full_cov_stringent_max_allele_w_coverage.results.txt", header=T, sep="\t", stringsAsFactors = F)
bad_sample_list<-subset(dat.dbscan.counts, outlier_count>148)$sample

#REad in dbscan output
dat.dbscan<-read.delim("/dbscan/eh_v5.dbscan_eps2.results_w_alleles_full_cov_stringent_max_allele_w_coverage.parsed.txt", header=T, sep="\t", stringsAsFactors = F)
dat.dbscan$str_id<-paste0(dat.dbscan$chr, "_", dat.dbscan$pos_start, "_", dat.dbscan$pos_end)
dat.dbscan<-subset(dat.dbscan, !sample%in%bad_sample_list)
dat.dbscan$cohort<-str_split_fixed(dat.dbscan$sample, "\\.", Inf)[,2]
dat.dbscan<-subset(dat.dbscan, cohort%in%c("ADC", "ROS", "MAP"))

#Generate genomic file of STR expansions seen in AD but not control samples; call this the foreground
case_str_list<-subset(dat.dbscan, AD==1)$str_id
control_str_list<-subset(dat.dbscan, AD==0)$str_id
case_str_list<-case_str_list[!case_str_list%in%control_str_list]
dat.counts<-data.frame(table(case_str_list), stringsAsFactors = F)
case_recurrent_str_list<-subset(dat.counts, Freq>=1)$case_str_list
dat.foreground<-data.frame(chr=NA, pos_start=NA, pos_end=NA, str_id=case_recurrent_str_list, stringsAsFactors=F)
dat.foreground$chr<-str_split_fixed(case_recurrent_str_list, "\\_", 3)[,1]
dat.foreground$pos_start<-as.numeric(str_split_fixed(case_recurrent_str_list, "\\_", 3)[,2])
dat.foreground$pos_end<-as.numeric(str_split_fixed(case_recurrent_str_list, "\\_", 3)[,3])
dat.foreground.se<-makeGRangesFromDataFrame(dat.foreground, seqnames = "chr", start.field = "pos_start", end.field = "pos_end",keep.extra.columns=TRUE) #Make into GRanges format

#Generate list of all STRs tested as background
dat.background<-fread("eh.v5_w_gangstr.v13.polymorphic.bed", sep="\t", header=F, data.table=F)
#dat.background<-fread(args[2], sep="\t", header=F, data.table=F)
dat.background$str_id<-paste0(dat.background$V1, "_", dat.background$V2, "_", dat.background$V3)
dat.background.se<-makeGRangesFromDataFrame(dat.background, seqnames = "V1", start.field = "V2", end.field = "V3",keep.extra.columns=TRUE) #Make into GRanges format


#Subset out segmental duplication regions from background
dat.segdup<-fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz", data.table=F, header=F, stringsAsFactors = F)
dat.segdup<-dat.segdup[,c(2:4)]
names(dat.segdup)<-c("chr", "start", "stop")
dat.segdup.gr<-makeGRangesFromDataFrame(dat.segdup, seqnames.field = "chr", start.field = "start", end.field = "stop")
overlap_index<-data.frame(findOverlaps(dat.background.se, dat.segdup.gr, type="any"))$queryHits
dat.background$index<-seq(1,nrow(dat.background))
dat.background$seg_dup<-ifelse(dat.background$index%in%overlap_index, 1, 0)
dat.background<-subset(dat.background, seg_dup==0, select=-c(seg_dup,index))
dat.background.se<-makeGRangesFromDataFrame(dat.background,
                                            seqnames = names(dat.background)[1], start.field = names(dat.background)[2], end.field = names(dat.background)[3],
                                            keep.extra.columns=TRUE) #Make into GRanges format




#run TE Fisher's exact test
dat.all<-data.frame(trial=numeric(), foreground_overlap=numeric(), background_overlap=numeric(), 
                    foreground_prop=numeric(), background_prop=numeric(), 
                    te=character(), dist=numeric(), stringsAsFactors = F)

te_list<-c("SVA","Alu","LINE2","LINE1","HERV")

for(s in te_list){

    #Read in repeat masker file and subset for TE class of interest
    dat.annot<-read.delim(paste0("/project/jcreminslab/guomic_projects/ad_str/eh_0822/enrichment/annotations/repeat_masker.",s, "_", d, "bp.hg38.bed"), header=F, sep="\t", stringsAsFactors = F)
    dat.annot$str_id<-paste0(dat.annot$V1, "_", dat.annot$V2, "_", dat.annot$V3)
    dat.annot.se<-makeGRangesFromDataFrame(dat.annot, seqnames = "V1", start.field = "V2", end.field = "V3",keep.extra.columns=TRUE) #Make into GRanges format
      
    #Run permutation testing
    n<-0
    
    dat.summary<-data.frame(trial=seq(0,n,1), foreground_overlap=0, background_overlap=0,stringsAsFactors = F)
    for(i in 1:nrow(dat.summary)){
      if(dat.summary[i,]$trial==0){
        dat.random.se<-dat.annot.se
      }else{
       # dat.random.se<-randomizeRegions(dat.annot.se, genome="hg38", per.chromsome=T)
        dat.random.se<-randomizeRegions(dat.annot.se, genome="hg38", per.chromsome=F, mask=dat.segdup.gr)
      }
      dat.summary[i,]$foreground_overlap<-length(unique(data.frame(subsetByOverlaps(dat.foreground.se, dat.random.se, type="any"), stringsAsFactors = F)$str_id)) #
      dat.summary[i,]$background_overlap<-length(unique(data.frame(subsetByOverlaps(dat.background.se, dat.random.se, type="any"), stringsAsFactors = F)$str_id))
    }
    dat.summary$foreground_prop<-dat.summary$foreground_overlap/nrow(dat.foreground)
    dat.summary$background_prop<-dat.summary$background_overlap/nrow(dat.background)
    dat.summary$te<-s
    dat.summary$dist<-d
    dat.all<-rbind(dat.all, dat.summary)
  }
}


#run ENCODE tissue enrichment for hippocampus
dat.all<-data.frame(trial=numeric(), foreground_overlap=numeric(), background_overlap=numeric(), 
                    foreground_prop=numeric(), background_prop=numeric(), 
                    tissue=character(),stringsAsFactors = F)
#encode_list<-"ENCFF886QLM"
encode_list<-c("h3k27ac","h3k27me3","h3k36me3","h3k4me1","h3k4me3","h3k9ac","h3k9me3", "h3k4me3_h3k27ac")
for(s in encode_list){
  dat.annot<-read.delim(paste0("/project/jcreminslab/guomic_projects/ad_str/eh_0822/annot/encode/", s, "/hippocampus.", s, ".encode.merged.bed"), header=F, sep="\t", stringsAsFactors = F)
  
    dat.annot$str_id<-paste0(dat.annot$V1, "_", dat.annot$V2, "_", dat.annot$V3)
    dat.annot.se<-makeGRangesFromDataFrame(dat.annot, seqnames = "V1", start.field = "V2", end.field = "V3",keep.extra.columns=TRUE) #Make into GRanges format
    

    #Run permutation testing
    n<-0
    
    dat.summary<-data.frame(trial=seq(0,n,1), foreground_overlap=0, background_overlap=0,stringsAsFactors = F)
    for(i in 1:nrow(dat.summary)){
      if(dat.summary[i,]$trial==0){
        dat.random.se<-dat.annot.se
      }else{
        #dat.random.se<-randomizeRegions(dat.annot.se, genome="hg38", per.chromsome=F)
        dat.random.se<-randomizeRegions(dat.annot.se, genome="hg38", per.chromsome=F, mask=dat.segdup.gr)
      }
      dat.summary[i,]$foreground_overlap<-length(unique(data.frame(subsetByOverlaps(dat.foreground.se, dat.random.se, type="any"), stringsAsFactors = F)$str_id)) #
      dat.summary[i,]$background_overlap<-length(unique(data.frame(subsetByOverlaps(dat.background.se, dat.random.se, type="any"), stringsAsFactors = F)$str_id))
    }
    dat.summary$foreground_prop<-dat.summary$foreground_overlap/nrow(dat.foreground)
    dat.summary$background_prop<-dat.summary$background_overlap/nrow(dat.background)
    dat.summary$tissue<-s
    dat.all<-rbind(dat.all, dat.summary)
}
write.table(dat.all, "eh_v5.dbscan_eps2.cases.E071_histone_enrichment.txt", row.names = F, col.names = T, sep="\t", quote=F)



#run H3K9me3 tissue enrichment
dat.all<-data.frame(trial=numeric(), foreground_overlap=numeric(), background_overlap=numeric(), 
                    foreground_prop=numeric(), background_prop=numeric(), 
                    encode_id=character(),stringsAsFactors = F)
#encode_list<-list.files("/project/jcreminslab/guomic_projects/ad_str/encode/h3k9me3/test", pattern="p0.9999.merge20kb_domain100kb.bed")
encode_list<-list.files("/project/jcreminslab/guomic_projects/ad_str/encode/h3k9me3/test", pattern="hippocampus")
encode_list<-encode_list[grepl(".rseg", encode_list)]
encode_list<-encode_list[grepl(".bed", encode_list)]

for(s in encode_list){
  dat.annot<-read.delim(paste0("/project/jcreminslab/guomic_projects/ad_str/encode/h3k9me3/test/",s), header=F, sep="\t", stringsAsFactors = F)
#  dat.annot<-read.delim(paste0("/project/jcreminslab/guomic_projects/ad_str/eh_0822/annot/encode/",s), header=F, sep="\t", stringsAsFactors = F)
  
  
  dat.annot$str_id<-paste0(dat.annot$V1, "_", dat.annot$V2, "_", dat.annot$V3)
  dat.annot.se<-makeGRangesFromDataFrame(dat.annot, seqnames = "V1", start.field = "V2", end.field = "V3",keep.extra.columns=TRUE) #Make into GRanges format
  
  #Run permutation testing
  n<-0
  dat.summary<-data.frame(trial=seq(0,n,1), foreground_overlap=0, background_overlap=0,stringsAsFactors = F)
  for(i in 1:nrow(dat.summary)){
    if(dat.summary[i,]$trial==0){
      dat.random.se<-dat.annot.se
    }else{
      dat.random.se<-randomizeRegions(dat.annot.se, genome="hg38", per.chromsome=F, mask=dat.segdup.gr)
    }
    dat.summary[i,]$foreground_overlap<-length(unique(data.frame(subsetByOverlaps(dat.foreground.se, dat.random.se, type="any"), stringsAsFactors = F)$str_id)) #
    dat.summary[i,]$background_overlap<-length(unique(data.frame(subsetByOverlaps(dat.background.se, dat.random.se, type="any"), stringsAsFactors = F)$str_id))
  }
  dat.summary$foreground_prop<-dat.summary$foreground_overlap/nrow(dat.foreground)
  dat.summary$background_prop<-dat.summary$background_overlap/nrow(dat.background)
  dat.summary$encode_id<-s
  dat.all<-rbind(dat.all, dat.summary)
}


#Calculate p-values using Fisher's exact test
dat.all$pval<-1
dat.all$or<-1
for(i in 1:nrow(dat.all)){
  a<-dat.all[i,]$foreground_overlap
  b<-round(dat.all[i,]$foreground_overlap/dat.all[i,]$foreground_prop)-a
  c<-dat.all[i,]$background_overlap
  d<-round(dat.all[i,]$background_overlap/dat.all[i,]$background_prop)-c
  if(a>0 & b>0 & c>0 & d>0){
  dat.all[i,]$pval<-fisher.test(rbind(c(a,c), c(b,d)), alternative="greater")$p
  dat.all[i,]$or<-as.numeric(fisher.test(rbind(c(a,c), c(b,d)), alternative="greater")$estimate)
  }
}
dat.all<-dat.all[order(dat.all$pval),]
write.table(dat.all, "eh_v5.dbscan_eps2.cases_recurrent1.all_tissues_chromhmm_enrichment.txt", row.names = F, col.names = T, sep="\t", quote=F)
