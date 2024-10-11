library(stringr)

sample_list<-list.files("/vcf/")

for(s in sample_list){
  dat.bed<-read.delim("eh.v5_w_gangstr.v13.polymorphic.bed", header=F, sep="\t", stringsAsFactors = F)
  names(dat.bed)<-c("chr", "pos_start", "pos_end")
  dat.bed<-subset(dat.bed, chr!="chrX")
  dat.bed$index<-seq(1, nrow(dat.bed))
  
  dat.vcf<-data.frame(chr=character(), pos_start=numeric(), gt=numeric(), stringsAsFactors = F)
  for(c in c(1:22)){
    
    #Read in vcf and obtain genotypes
    filename<-paste0("/vcf/", s, "/", s, ".eh-v5.chr", c, ".vcf.gz")
    dat.temp<-read.delim(filename, comment.char = "#", header=F, sep="\t")
    names(dat.temp)[c(1:2)]<-c("chr", "pos_start")
    dat.temp$gt<-str_split_fixed(dat.temp$V10, "\\:", 8)[,3]
    dat.temp<-subset(dat.temp, select=c(chr, pos_start, gt))
    
    #For STRs with multiple entries, take the one with the longest allele
    dat.temp$long<-as.numeric(str_split_fixed(dat.temp$gt, "\\/", 2)[,2])
    dat.temp<-dat.temp[order(dat.temp$pos_start, -dat.temp$long),]
    dat.temp<-dat.temp[!duplicated(dat.temp$pos_start),]
    dat.temp<-subset(dat.temp, select=-c(long))
    dat.vcf<-rbind(dat.vcf, dat.temp)
  }

  #Merge genotypes with bed file and reorder
  dat.bed<-merge(dat.bed, dat.vcf, by=c("chr", "pos_start"), all.x=T, all.y=F)
  dat.bed<-dat.bed[order(dat.bed$index),]
  
  #output data
  dat.bed<-subset(dat.bed, select=c(gt))
  names(dat.bed)<-s
  write.table(dat.bed, paste0( s, ".eh-v5.genotype.txt"), row.names = F, col.names = T, sep="\t", quote=F)
}
