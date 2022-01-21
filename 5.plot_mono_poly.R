library(RColorBrewer)
library(ggplot2)
library(reshape2)

projectdir="/share/ScratchGeneral/nenbar/projects/canfam"
scriptsdir=paste0(projectdir,"/scripts")
resultsdir=paste0(projectdir,"/results")
figuredir=paste0(projectdir,"/figures")
system(paste0("mkdir -p ",figuredir))

populations<-read.table(paste0(projectdir,"/annotation/popmap_0702CR30"),sep="\t")

pops<-gsub("Galaxias ","",populations$V2)
popsT<-as.data.frame(table(pops))


results<-list()
missing<-c("0.8","0.9","1")
for(missingness in missing){

       #the output file consists of two joint tables with different number of rows
       #1. extract table 1 for polymorphic loci and find the number of samples
       dataPoly<-read.table(paste0("../results/",missingness,"_miss/populations.sumstats_summary.tsv.backup"),sep="\t",header=T,fill=T,comment.char = "",skip=1)
       #find the first comment line
       commentN<-min(grep("#",dataPoly[,1]))
       nSamples<-commentN-1
       dataPoly<-dataPoly[1:nSamples,]
       colnames(dataPoly)[1]<-"Pop.ID"
       #subselect only certain columns
       dataPolyS<-dataPoly[,c("Pop.ID","Obs_Het","Exp_Het")]
       #2. exctract table 2 for monorphic loci only
       dataMono<-read.table(paste0("../results/",missingness,"_miss/populations.sumstats_summary.tsv.backup"),sep="\t",header=T,fill=T,comment.char = "",skip=commentN+2)
       colnames(dataMono)[1]<-"Pop.ID"
       dataMonoS<-dataMono[,c("Pop.ID","Obs_Het","Exp_Het")]

       results[[missingness]]<-rbind(dataPolyS,dataMonoS)
}

df<-do.call("rbind",results)
dataM<-melt(df,id.vars="Pop.ID")
dataM$missingness<-rep(missing,each=nSamples*4)
dataM$type<-rep(c("poly","mono"),each=nSamples*2)

#cols<-c(brewer.pal(3,"Reds")[1:2],brewer.pal(3,"Blues")[1:2],brewer.pal(3,"Greens")[1:2])

dataM$Pop.ID<-gsub("Galaxias ","",dataM$Pop.ID)

dfM<-merge(dataM,popsT,by.x="Pop.ID",by.y="pops")
dfM$Pop.ID<-paste0(dfM$Pop.ID,"\n(",dfM$Freq,")")
colnames(dfM)<-c("population","class","heterozygosity","missingness","type","pop.size")
dfM$class<-gsub("_.*,","",dfM$class)
dfM$heterozygosity<-as.numeric(dfM$heterozygosity)

pdf(paste0(figuredir,"/missingness_mono_poly.pdf"),width=12,height=8)
p<-ggplot(dfM,aes(x=population,y=heterozygosity,color=class,shape=type))
p<-p+geom_point()
p<-p+theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
p<-p+scale_y_log10()
p
dev.off()
