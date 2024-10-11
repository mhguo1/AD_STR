#This code performs a hypergeometric Fisher's exact test for number of STR expansion carriers at different tract length thresholds

library(data.table)
library(stringr)

dat.annot<-read.delim("eh.v5_w_gangstr.v13.polymorphic.bed", header=T, sep="\t", stringsAsFactors = F) #read in bed file of STRs

#Read in merged STR genotypes
dat<-fread("EH_v5.merged.max_allele.txt.gz", header=T, sep="\t", data.table=F)
dat<-cbind(dat.annot, dat[,4:ncol(dat)])

#read in coverage file
dat.site_cov<-fread("EH_v5.merged.coverage.txt.gz", header=T, sep="\t", data.table=F)
dat.site_cov<-cbind(dat.annot, dat.site_cov[,4:ncol(dat.site_cov)])

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

#subset on samples of interest, coverage, and PC coordinates
dat.pheno<-dat.pheno[grepl("ADC", dat.pheno$SUBJID),]
dat.pheno<-subset(dat.pheno,PC1< (-0.0037) & PC2<0.02 & COV_avg>30 & COV_avg<50 & Platform=="HiSeqX")

#Generate case and control columns
dat_names<-names(dat)
dat_names<-str_split_fixed(dat_names, ".B", 2)[,1]

case_cols<-which(dat_names%in%dat.pheno[dat.pheno$AD==1,]$SUBJID)
control_cols<-which(dat_names%in%dat.pheno[dat.pheno$AD==0,]$SUBJID)
last_col<-ncol(dat)

##Remove STRs with no standard deviation
dat$sd<-apply(dat[,c(6:last_col)], 1, sd, na.rm=T)
dat<-subset(dat, sd>0, select=-c(sd))

#Remove samples with genotyping rate < 90%
geno_rate<-function(temp_vector){
  length(temp_vector[!is.na(temp_vector)])/length(temp_vector)
}
dat$case_geno_rate<-apply(dat[,case_cols], 1, geno_rate)
dat$control_geno_rate<-apply(dat[,control_cols], 1, geno_rate)
dat<-subset(dat, case_geno_rate>0.9 & control_geno_rate>0.9, select=-c(case_geno_rate, control_geno_rate))

#Calculate change in tract length
for(c in 6:ncol(dat)){
  dat[,c]<-as.numeric(dat[,c])-dat$ref
}

#Create dat.test output file
dat.test<-dat[,c(1:5)]
dat.test$index<-seq(1,nrow(dat.test)) #Create row index for faster testing later

#calculate number of individuals with non NA STR genotypes
countif<-function(temp_vector){
  length(temp_vector[!is.na(temp_vector)])
}
dat.test$case_count<-apply(dat[,case_cols], 1, countif)
dat.test$control_count<-apply(dat[,control_cols], 1, countif)

#Perform Fisher's exact test at different STR tract length thresholds (parameter s)
for(s in c(1,5,10,20)){
  #s is a cutoff for expansion length
  
  #Generate case and control counts for each expansion
  countif<-function(temp_vector){
    length(temp_vector[temp_vector>=s & !is.na(temp_vector)])
  }
  dat.test$case_exp_temp<-apply(dat[,case_cols], 1, countif) #Calculate number of case individuals with STR tract length exceeding threshold
  dat.test$control_exp_temp<-apply(dat[,control_cols], 1, countif) #Calculate number of control  individuals with STR tract length exceeding threshold
  dat.test$pval_temp<-1
  
  #Calculate Fisher's exact test p-value
  nonzero_cols<-subset(dat.test, case_exp_temp>0 | control_exp_temp>0)$index
  for(i in nonzero_cols){
    dat.test[i,]$pval_temp<-fisher.test(rbind(c(dat.test[i,]$case_exp_temp,dat.test[i,]$case_count-dat.test[i,]$case_exp_temp),
                                              c(dat.test[i,]$control_exp_temp,dat.test[i,]$control_count-dat.test[i,]$control_exp_temp)))$p.value
  }
  
  names(dat.test)<-gsub("case_exp_temp", paste0("case_exp", s), names(dat.test))
  names(dat.test)<-gsub("control_exp_temp", paste0("control_exp", s), names(dat.test))
  names(dat.test)<-gsub("pval_temp", paste0("pval_exp", s), names(dat.test))
}

write.table(dat.test, "eh_v5.polymorphic.fisher.adc_rosmap.results.txt", row.names=F, col.names = T, sep="\t", quote=F)
