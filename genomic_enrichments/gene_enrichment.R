library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(stringr)
library(biomaRt)
library(ggplot2)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#Read in DBSCAN results
dat<-read.delim("eh_v5.dbscan_eps2.results_w_alleles_full_cov_stringent_max_allele_w_coverage.hiseqx.gencode_annotation.txt", 
                header=T, sep="\t", stringsAsFactors = F)


###Perform enrichment analyses for genes within window
#Find genes within 250kb window
case_gene_list<-unlist(strsplit(paste(subset(dat, case_expansion_count>0 & abs(tss_distance)<250000 &  !is.na(tss_gene))$tss_gene, collapse=","), "\\,"))
background_gene_list<-unique(unlist(strsplit(paste(subset(dat, abs(tss_distance)<250000 & !is.na(tss_gene))$tss_gene, collapse=","), "\\,")))

foreground_entrez <- mapIds(org.Hs.eg.db, keys = case_gene_list, column = "ENTREZID", keytype = "SYMBOL")
background_entrez <- mapIds(org.Hs.eg.db, keys = background_gene_list, column = "ENTREZID", keytype = "SYMBOL")

# Perform enrichment analysis
enrich_result <- enrichGO(gene         = foreground_entrez,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "ENTREZID",
                          ont          = "ALL",  # Specify ontology ("BP", "MF", or "CC")
                          universe     = background_entrez,
                          pvalueCutoff = 0.05,   # Adjust as needed
                          qvalueCutoff = 1.0)    # Adjust as needed
dat.results<-data.frame(enrich_result, stringsAsFactors = F)
dat.results<-subset(dat.results, select=-c(geneID))

#Plot enrichment of genes within 250kb window (Figure 5A)
dat.plot<-dat.results[order(dat.results$pvalue),]
dat.plot$Description<-factor(dat.plot$Description,levels=rev(unique(dat.plot$Description)))
dat.plot$or<-(as.numeric(str_split_fixed(dat.plot$GeneRatio, "\\/", 2)[,1])/as.numeric(str_split_fixed(dat.plot$GeneRatio, "\\/", 2)[,2]))/(as.numeric(str_split_fixed(dat.plot$BgRatio, "\\/", 2)[,1])/as.numeric(str_split_fixed(dat.plot$BgRatio, "\\/", 2)[,2]))

ggplot(dat.plot[c(1:10),],aes(x=Description,y=-log10(pvalue)))+ 
  geom_point(cex=4)+coord_flip(ylim=c(5,13))+
  theme_classic()+scale_color_brewer(palette = "Set2") #+scale_size(breaks=c(1,2))



###Perform enrichment analyses for STRs within gene bodies
#Find genes with STRs within gene bodies
case_gene_list<-unique(unlist(strsplit(paste(subset(dat, case_expansion_count>0 & control_expansion_count==0 &!is.na(gene))$gene, collapse=","), "\\,")))
background_gene_list<-unique(unlist(strsplit(paste(subset(dat,  !is.na(gene))$gene, collapse=","), "\\,")))

foreground_entrez <- mapIds(org.Hs.eg.db, keys = case_gene_list, column = "ENTREZID", keytype = "SYMBOL")
background_entrez <- mapIds(org.Hs.eg.db, keys = background_gene_list, column = "ENTREZID", keytype = "SYMBOL")

# Perform enrichment analysis
enrich_result <- enrichGO(gene         = foreground_entrez,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "ENTREZID",
                          ont          = "ALL",  # Specify ontology ("BP", "MF", or "CC")
                          universe     = background_entrez,
                          pvalueCutoff = 0.05,   # Adjust as needed
                          qvalueCutoff = 1.0)    # Adjust as needed
dat.results<-data.frame(enrich_result, stringsAsFactors = F)
dat.results<-subset(dat.results, select=-c(geneID))

#Plot gene body enrichments (Figure 5B)
dat.plot<-dat.results[order(dat.results$pvalue),]
dat.plot$Description<-factor(dat.plot$Description,levels=rev(unique(dat.plot$Description)))
dat.plot$or<-(as.numeric(str_split_fixed(dat.plot$GeneRatio, "\\/", 2)[,1])/as.numeric(str_split_fixed(dat.plot$GeneRatio, "\\/", 2)[,2]))/(as.numeric(str_split_fixed(dat.plot$BgRatio, "\\/", 2)[,1])/as.numeric(str_split_fixed(dat.plot$BgRatio, "\\/", 2)[,2]))

ggplot(dat.plot[c(1:10),],aes(x=Description,y=-log10(pvalue)))+ 
  geom_point(cex=4)+coord_flip(ylim=c(5,13))+
  theme_classic()+scale_color_brewer(palette = "Set2")
