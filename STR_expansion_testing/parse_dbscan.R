#Parse dbscan output
library(stringr)
library(data.table)
library(GenomicRanges)

#Parse output from DBSCAN script
dat.dbscan<-read.delim(paste0("eh_v5.dbscan_eps2.results_w_alleles_full_cov_stringent_max_allele_w_coverage.txt"), stringsAsFactors = F, sep="\t", header=T)
dat.dbscan<-subset(dat.dbscan, !is.na(outliers))

#Subset out segmental duplication regions#
dat.dbscan$str_id<-paste0(dat.dbscan$chr, "_", dat.dbscan$pos_start, "_", dat.dbscan$pos_end)
dat.dbscan.gr<-makeGRangesFromDataFrame(dat.dbscan, seqnames.field = "chr", start.field = "pos_start", end.field = "pos_end", keep.extra.columns = T)

dat.segdup<-fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz", data.table=F, header=F, stringsAsFactors = F)
dat.segdup<-dat.segdup[,c(2:4)]
names(dat.segdup)<-c("chr", "start", "stop")
dat.segdup.gr<-makeGRangesFromDataFrame(dat.segdup, seqnames.field = "chr", start.field = "start", end.field = "stop")

overlap_index<-data.frame(findOverlaps(dat.dbscan.gr, dat.segdup.gr, type="any"))$queryHits
dat.dbscan$index<-seq(1,nrow(dat.dbscan))
dat.dbscan$seg_dup<-ifelse(dat.dbscan$index%in%overlap_index, 1, 0)
seg_dups<-subset(dat.dbscan, seg_dup==1)$str_id

dat.dbscan<-subset(dat.dbscan, seg_dup==0)

#Summaryize dbscan output
dat.out<-data.frame(chr=character(), pos_start=numeric(), pos_end=numeric(), motif=character(), sample=character(), outlier_length=numeric(), stringsAsFactors = F)
for(i in c(1:nrow(dat.dbscan))){
  dat.temp<-data.frame(chr=dat.dbscan[i,]$chr, pos_start=dat.dbscan[i,]$pos_start, pos_end=dat.dbscan[i,]$pos_end, motif=dat.dbscan[i,]$ru, ref=dat.dbscan[i,]$ref,
                       sample=as.vector(str_split_fixed(dat.dbscan[i,]$outliers, "\\;", Inf)),outlier_length=as.vector(as.numeric(str_split_fixed(dat.dbscan[i,]$outlier_size, "\\;", Inf))), 
                       stringsAsFactors = F)
  dat.out<-rbind(dat.temp, dat.out)
}

#Generate outlier counts per sample
dat.out.sample<-data.frame(table(dat.out$sample), stringsAsFactors = F)
names(dat.out.sample)<-c("sample", "outlier_count")

#Read in PCA file (downloaded from NIAGADS)
dat.pc<-read.delim("/pheno/17k_phenotype_pc_der_20220209.fix_PCR.csv", header=T, sep=",", stringsAsFactors = F)
dat.pc<-subset(dat.pc, !is.na(AD) & !is.na(age) & PCR_Free=="Yes" & Eth=="NonHispanicWhite", select=c(ID, PC1, PC2, PC3, sex, age, AD, Sequencing_Center, Platform))

#Read in sample coverage file (downloaded from NIAGADS)
dat.cov<-read.delim("/pheno/17k_coverage.txt", header=T, sep="\t", stringsAsFactors = F)
dat.cov<-subset(dat.cov, select=c(Sample, COV_avg))
names(dat.cov)[1]<-"ID"
dat.pc<-merge(dat.pc, dat.cov, by="ID", all.x=F, all.y=F)
dat.pc$SUBJID<-paste(str_split_fixed(dat.pc$ID, "\\-", Inf)[,1],str_split_fixed(dat.pc$ID, "\\-", Inf)[,2],str_split_fixed(dat.pc$ID, "\\-", Inf)[,3], sep=".")
dat.pc<-subset(dat.pc, select=-c(ID))

#read in phenotype file (downloaded from NIAGADS)
dat.pheno<-read.delim("ADSPCaseControlPhenotypes_DS_2021.02.19.v2_ALL.txt")
dat.pheno$apoe<-as.numeric(substr(dat.pheno$APOE, 1,1))+as.numeric(substr(dat.pheno$APOE, 2,2))
dat.pheno<-subset(dat.pheno, select=c(SUBJID, apoe))
dat.pheno$SUBJID<-gsub("\\-", "\\.", dat.pheno$SUBJID)

#Merge and process phenotype files
dat.pheno<-merge(dat.pc, dat.pheno, by="SUBJID", all.x=F, all.y=F)
dat.pheno$AD<-factor(dat.pheno$AD, levels=c(0,1))
dat.pheno$Sequencing_Center<-factor(dat.pheno$Sequencing_Center)
dat.pheno$Platform<-factor(dat.pheno$Platform)

dat.pheno<-subset(dat.pheno,PC1< (-0.0037) & PC2<0.02 & COV_avg>30 & COV_avg<50 & Platform=="HiSeqX")
names(dat.pheno)[1]<-"sample"

#Add back in samples with no expansions identified from dbscan
sample_list<-read.delim("all_sample.list.txt", header=F, sep="\t", stringsAsFactors = F)$V1
if(length(sample_list[!sample_list%in%dat.out.sample$sample])>0){
  dat.temp<-data.frame(sample=sample_list[!sample_list%in%dat.out.sample$sample], outlier_count=0, stringsAsFactors = F)
  dat.out.sample<-rbind(dat.out.sample, dat.temp)
}

#Write out each expansion identified in DBSCAN along with sample information
dat.out<-merge(dat.out, dat.pheno, all.x=T, all.y=F, by="sample")
write.table(dat.out, "eh_v5.dbscan_eps2.results_w_alleles_full_cov_stringent_max_allele_w_coverage.parsed.txt", row.names=F, col.names=T, sep="\t", quote=F)

#Write out counts of number of expansions for each sample
dat.out.summary<-merge(dat.out.sample, dat.pheno, all.x=T, all.y=F, by="sample")
dat.out.summary$cohort<-str_split_fixed(dat.out.summary$sample, "\\.", Inf)[,2]
write.table(dat.out.summary, "eh_v5.dbscan_eps2.results_w_alleles_full_cov_stringent_max_allele_w_coverage.results.txt", row.names=F, col.names=T, sep="\t", quote=F)
