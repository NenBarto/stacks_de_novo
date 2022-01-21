
projectdir="/share/ScratchGeneral/nenbar/projects/canfam"
scriptsdir=paste0(projectdir,"/scripts")
annotationdir=paste0(projectdir,"/annotation/samples")

files<-list.files(annotationdir,pattern="targets",full.names=T)
results<-list()
for(file in files){
	results[[file]]<-read.csv(file)
}
targets<-do.call("rbind",results)

data<-read.csv(paste0(annotationdir,"/Final_Dingo_Data_4Nov2021.csv"))


#clean targets
targets$id<-paste0("a",targets$targetid)
targets$cleanIDs<-targets$genotype
targets$cleanIDs<-gsub("_rep.*","",targets$cleanIDs)
merged<-merge(targets,data,by.x="cleanIDs",by.y="Sample.ID")
write.table(merged,"combined_annotation.txt",quote=F,sep="\t")

popmap<-merged[,c("targetid","State")]
popmap$State[merged$LACP.Region=="Mallee"]<-"Mallee"
popmap$State[popmap$State=="Qld"]<-"QLD"
write.table(popmap,file=paste0(paste0(projectdir,"/annotation/popmap")),quote=F,row.names=F,col.names = F,sep="\t")