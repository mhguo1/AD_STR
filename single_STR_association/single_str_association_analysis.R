###Run single STR association analysis
library(data.table)
library(stringr)
library(RNOmni)

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

#Subset on samples shared between phenotype files and STR genotyping files
names(dat)[6:ncol(dat)]<-paste(str_split_fixed(names(dat)[6:ncol(dat)], "\\.", Inf)[,1], str_split_fixed(names(dat)[6:ncol(dat)], "\\.", Inf)[,2], str_split_fixed(names(dat)[6:ncol(dat)], "\\.", Inf)[,3], sep=".")
names(dat.site_cov)[6:ncol(dat.site_cov)]<-paste(str_split_fixed(names(dat.site_cov)[6:ncol(dat.site_cov)], "\\.", Inf)[,1], str_split_fixed(names(dat.site_cov)[6:ncol(dat.site_cov)], "\\.", Inf)[,2], str_split_fixed(names(dat.site_cov)[6:ncol(dat.site_cov)], "\\.", Inf)[,3], sep=".")
dat.pheno<-subset(dat.pheno, SUBJID%in%names(dat))

#Generate case and control sample lists
case_list<-subset(dat.pheno, AD==1 & SUBJID%in%names(dat))$SUBJID
control_list<-subset(dat.pheno, AD==0 & SUBJID%in%names(dat))$SUBJID
dat<-dat[,c(1:5, which(names(dat)%in%c(case_list, control_list)))]

#Process STR site coverage file
dat.site_cov<-dat.site_cov[,c(names(dat))]
dat.site_cov$str<-paste0(dat.site_cov[,c(1:3)], collapse="_")


#Helper function to calculate mode of STR tract lengths
Mode <- function(x){ 
  a = table(as.numeric(x)) # x is a vector
  return(as.numeric(names(a[which.max(a)])))
}

#Helper function to calculate number of non-reference STR alleles
non_ref_count<-function(temp_vector){ 
  temp_vector<-as.numeric(temp_vector)
  length(temp_vector[!is.na(temp_vector) & temp_vector!=Mode(temp_vector)])
}


#Generate output file
dat.test<-dat[,c(1:5)]
dat.test$effect<-0
dat.test$logistic_p<-1
dat.test$invnorm_p<-1
dat.test$mannwhitney_p<-1
dat.test$mean_control<-0
dat.test$mean_ad<-0
dat.test$sd<-0
dat.test$num_allele<-0
dat.test$alleles_ad<-NA
dat.test$alleles_control<-NA

#Run STR associations.
for(i in 1:nrow(dat.test)){
  if(i%%1000==0){print(i)}
    dat.temp<-data.frame(t(dat[i,c(6:ncol(dat))]), stringsAsFactors = F) #Generate temporary dataframe of genotypes for STR
    str<-paste0(dat[i,c(1:3)], collapse="_")
  
    names(dat.temp)[1]<-"str"
    non_ref<-non_ref_count(dat.temp$str) #Calculate number of non-reference STR alleles
    dat.temp$SUBJID<-rownames(dat.temp)
    
    dat.site_cov_temp<-subset(dat.site_cov, str==str, select=-c(str)) #Generate temporary dataframe of local coverage for STR
    dat.site_cov_temp<-data.frame(t(dat.site_cov_temp[i,c(6:ncol(dat.site_cov_temp))]), stringsAsFactors = F)
    dat.temp$site_cov<-dat.site_cov_temp[,1]  
    
    dat.temp<-merge(dat.temp, dat.pheno, by="SUBJID", all.x=F, all.y=F)
    dat.temp$str<-as.numeric(dat.temp$str)
    dat.temp<-subset(dat.temp, !is.na(str))  
    
    dat.temp$inv_norm<-RankNorm(as.numeric(dat.temp$str))   #Inverse normal transform STR genotypes
    dat.test[i,]$num_allele<-length(unique(dat.temp$str)) #Calculate number of unique alleles

    #Run association analysis
    if(non_ref>0 & nrow(dat.temp)>0){
      ad_table<-data.frame(table(subset(dat.temp, AD==1)$str))
      dat.test[i,]$alleles_ad<-paste(paste(ad_table[,1], ad_table[,2], sep=","), collapse="|") #Generate list of case STR alleles
      
      control_table<-data.frame(table(subset(dat.temp, AD==0)$str))
      dat.test[i,]$alleles_control<-paste(paste(control_table[,1], control_table[,2], sep=","), collapse="|") #Generate list of control STR alleles
      
      dat.test[i,]$mean_control<-mean(as.numeric(subset(dat.temp, AD==0)$str)) #Calculate mean STR tract length in controls
      dat.test[i,]$mean_ad<-mean(as.numeric(subset(dat.temp, AD==1)$str)) #Calculate mean STR tract length in cases
      dat.test[i,]$sd<-sd(as.numeric(dat.temp$str)) #Calculate standard deviation of STR tract lengths
      
      dat.test[i,]$logistic_p<-coef(summary(glm(AD ~ str+sex + PC1 + PC2 + PC3 + COV_avg + site_cov, data = dat.temp, family = "binomial")))[2,4] #Calculate p-value under logistic model using untransformed STR lengths
      dat.test[i,]$invnorm_p<-coef(summary(glm(AD ~ inv_norm+sex + PC1 + PC2 + PC3 + COV_avg + site_cov, data = dat.temp, family = "binomial")))[2,4] #Calculate p-value under logistic mdoel using inverse normal transformed STR lengths
      dat.test[i,]$mannwhitney_p<-wilcox.test(str~ AD, data = dat.temp)$p.value #Calculate p-value from non-parametric test
      dat.test[i,]$effect<-coef(summary(glm(AD ~ str+sex + PC1 + PC2 + PC3 +COV_avg + site_cov, data = dat.temp, family = "binomial")))[2,1] #Calculate effect size and direction under logistic model
    }
  }
write.table(dat.test, paste0("eh_v5.inv_norm.results_w_alleles_full_cov_stringent_max_allele_w_coverage.adc_rosmap.txt"), row.names=F, col.names = T, sep="\t", quote=F) #Write output file
