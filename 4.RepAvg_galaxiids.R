library(vcfR)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(ggdendro)
library(ggplot2)
library(dendextend)
library(dynamicTreeCut)

projectdir="/share/ScratchGeneral/nenbar/projects/canfam"
scriptsdir=paste0(projectdir,"/scripts")
annotationdir=paste0(projectdir,"/annotation/samples")
tabledir=paste0(projectdir,"/tables")
system(paste0("mkdir -p ",tabledir))

#import vcf file
missing<-c("1")
for(missingness in missing){
       #vcf.fn<-"../results/0702CR30.opt.stacksoptimised3_3/populations.snps.vcf.gz"
       vcf.fn<-paste0("../results/",missing,"_miss/populations.snps.vcf.gz")
       snpgdsVCF2GDS(vcf.fn, paste0("ccm_",missingness,".gds"))
       genofile <- snpgdsOpen(paste0("ccm_",missingness,".gds"))

       #get genotypes and sample IDs
       g <- read.gdsn(index.gdsn(genofile, "genotype"))
       sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

       #load in the annotation to get replicate data from "target" files from dartseq
       files<-list.files(annotationdir,pattern="targets_galaxiids",full.names=T)
       results<-list()
       for(file in files){
       	results[[file]]<-read.csv(file)
       }
       targets<-do.call("rbind",results)
       data<-read.csv(paste0(annotationdir,"/galaxiids_20210322_AWedited3Apr2021.csv"))

       #clean targets
       targets$id<-paste0("a",targets$targetid)
       targets$cleanIDs<-targets$genotype
       targets$cleanIDs<-gsub("_rep.*","",targets$cleanIDs)
       merged<-merge(targets,data,by.x="cleanIDs",by.y="Genotype")
       write.table(merged,paste0(tabledir,"/combined_annotation_0702CR30_",missingness,".txt"),sep="\t")

       #clean genotypes
       gS<-as.data.frame(g)
       gS$id<-sample.id
       gS<-gS[sample.id %in% as.character(merged$targetid),]
       mergeG<-merge(gS,merged[,c("cleanIDs","targetid")],by.x="id",by.y="targetid")

       #find the proportion of technical pairs for which the marker score is consistent
       #turn to matrix for faster processing
       mergeGMat<-as.matrix(mergeG[,-c(1,dim(mergeG)[2])])
       ids<-mergeG[,"cleanIDs"]

       #split matrix by ID
       mergeGL<-lapply(split(seq_along(ids), ids), #split indices by a
              function(m, ind) m[ind,], m = mergeGMat)[order(unique(ids))]

              #extract only duplicates

       #get only duplicate data
       isDuplicate<-sapply(mergeGL,function(x){class(x)[1]!="integer"})
       mergedDups<-mergeGL[isDuplicate]
       consistentCount<-sapply(mergedDups,function(x){apply(x,2,function(x){c(abs(max(x)-min(x)))==0})})

       #there were 80 duplicates, calculate sum of identical elements and divide by 80
       RepAvg<-rowSums(consistentCount)/dim(consistentCount)[2]
       snp.chromosome<-read.gdsn(index.gdsn(genofile, "snp.chromosome"))
       snp.position<-read.gdsn(index.gdsn(genofile, "snp.position"))
       df<-data.frame(chr=snp.chromosome,position=snp.position,RepAvg=RepAvg)
       write.table(df,paste0(tabledir,"/galaxiids_RepAvg_",missingness,".txt"),sep="\t",quote=F,row.names=F)
       #snpgdsClose(genofile)
}

ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2))
#ibs.hc$sample.id=paste0(gsub("Galaxias ","",data$V2),"_",data$V1)
data<-read.table("../annotation/popmap_0702CR30",sep="\t")
rv <- snpgdsCutTree(ibs.hc,samp.group=as.factor(gsub("Galaxias ","",data$V2)))

rvTemp=rv
#rvTemp$sample.id=paste0(gsub("Galaxias ","",data$V2),"_",data$V1)
#rvTemp$sample.id=paste0(gsub("Galaxias ","",data$V2))

pop<-factor(gsub("Galaxias ","",data$V2))

pdf("tree.pdf",width=20,height=8)
plot(rv$dendrogram, leaflab="none", main="Missingness 100%",xaxt="n",col="blue")
axis(2,cex.axis=1)
legend("topright", legend=levels(pop), col=rainbow(nlevels(pop)), pch=19, ncol=4)

dev.off()

hcl<-ibs.hc$hclust
hcl$labels<-gsub("Galaxias ","",data$V2)
#dhc <- as.dendrogram(hcl)
## Rectangular lines
#ddata <- dendro_data(dhc, type = "rectangle")
#
#dend<-as.dendrogram(hcl)
#ggd1 <- as.ggdend(dend)
#
pdf("tree.pdf",width=8,height=20)
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) +
  geom_text(data = label(ddata), 
              aes(x = x, y = y, label = label), hjust = -1, size = 1)
p
dev.off()

#df2<-data.frame(cluster=cutree(hcl,6),ids=factor(hcl$labels,levels=hcl$labels[hcl$order]))
#
#clusters<-data.frame(cluster=1:length(unique(pop)),pop=unique(pop))
#data$V2<-gsub("Galaxias ","",data$V2)
#dataM<-merge(data,clusters,by.x="V2",by.y="pop")
#df2M<-merge(df2,dataM,by.x="ids",by.y="V1")
#
#df2<-df2M[,c("cluster.y","ids")]
#colnames(df2)<-c("cluster","ids")
#
#p1<-ggdendrogram(hc, rotate=FALSE)+
#  theme(axis.title=element_blank(),
#        axis.ticks=element_blank(),
#        axis.text=element_blank(),
#        legend.position="none")
#
#
#p2<-ggplot(df2,aes(states,y=1,fill=factor(cluster)))+geom_tile()+
#  scale_y_continuous(expand=c(0,0))+
#  theme(axis.title=element_blank(),
#        axis.ticks=element_blank(),
#        axis.text=element_blank(),
#        legend.position="none")
#
#library(gridExtra)
#
#gp1<-ggplotGrob(p1)
#gp2<-ggplotGrob(p2)  
#
#maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
#gp1$widths[2:5] <- as.list(maxWidth)
#gp2$widths[2:5] <- as.list(maxWidth)
#
#pdf("tree.pdf",width=20,height=8)
#
#grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5))
#
#dev.off()#