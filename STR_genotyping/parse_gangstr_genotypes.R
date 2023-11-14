library(stringr)

args = commandArgs(trailingOnly=TRUE)

s=args[1] #path to vcf file
for(s in sample_list){
  dat.bed<-read.delim("gangstr.v13.polymorphic_w_eh.v5_offtarget.exp10_p05.top5_no_segdup.fast_2s.bed", header=F, sep="\t", stringsAsFactors = F)
  dat.bed<-dat.bed[,c(1:3)]
  names(dat.bed)<-c("chr", "pos_start", "pos_end")
 # dat.bed<-subset(dat.bed, chr!="chrX")
  dat.bed$index<-seq(1, nrow(dat.bed))
  

  #Read in vcf and obtain genotypes
  filename<-s
  dat.temp<-read.delim(filename, comment.char = "#", header=F, sep="\t")
  names(dat.temp)[c(1:2)]<-c("chr", "pos_start")
  dat.temp$gt<-str_split_fixed(dat.temp$V10, "\\:", Inf)[,4]
  dat.vcf<-subset(dat.temp, select=c(chr, pos_start, gt))
    
  #Merge genotypes with bed file and reorder
  dat.vcf<-merge(dat.bed, dat.vcf, by=c("chr", "pos_start"), all.x=T, all.y=F)
  dat.vcf<-dat.vcf[order(dat.vcf$index),]
  
  #output data
  dat.vcf<-subset(dat.vcf, select=c(gt))
  s<-gsub(".gangstr.v13.polymorphic_w_eh.v5_offtarget.genotype.vcf", "", s)
  names(dat.vcf)<-s
  write.table(dat.vcf, paste0( s, ".gangstr.v13.polymorphic_w_eh.v5_offtarget.genotype.txt"), row.names = F, col.names = T, sep="\t", quote=F)
}
