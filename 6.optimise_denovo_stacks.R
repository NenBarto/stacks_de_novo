
projectname="0702CR30"
projectdir="/share/ScratchGeneral/nenbar/projects/canfam"
scriptsdir=paste0(projectdir,"/scripts")
annotationdir=paste0(projectdir,"/annotation/samples")
rawdir=paste0(projectdir,"/raw/",projectname,"/combined")


#first extract the dart data
files<-list.files(rawdir,pattern="targets",full.names=T)
results<-list()
for(file in files){
	results[[file]]<-read.csv(file)
}
targets<-do.call("rbind",results)

#then get the population data
data<-read.csv(paste0(annotationdir,"/galaxiids_20210322_AWedited3Apr2021.csv"))


#clean targets
targets$id<-paste0("a",targets$targetid)
targets$cleanIDs<-targets$genotype
#targets$cleanIDs<-gsub("rep.*","",targets$cleanIDs)
merged<-merge(targets,data,by.x="cleanIDs",by.y="Genotype")
write.table(merged,paste0("combined_annotation_",projectname,".txt"),quote=F,sep="\t")

#create subpopulations
#eliminate small populations and duplicates
mergedS<-merged[!(merged$Species %in% c("Galaxias cann","Galaxias terenasus")),]
mergedS<-mergedS[!(grepl("rep",mergedS$genotype)),]
dupIDs<-table(mergedS$genotype)
dupIDs<-names(dupIDs[dupIDs>1])
mergedS<-mergedS[!(mergedS$genotype %in% dupIDs),]

mergedSL<-split(mergedS[,c("targetid","genotype","Species")],mergedS$Species)
set.seed(1)
mergedSLsubpop<-lapply(mergedSL,function(x){x[sample(1:dim(x)[1],5),]})
subpopdf<-do.call("rbind",mergedSLsubpop)

popmapShort<-subpopdf[,c("targetid","Species")]
popmapShort$Species="opt"
write.table(popmapShort,file=paste0(paste0(projectdir,"/annotation/popmapOptimise_",projectname)),quote=F,row.names=F,col.names = F,sep="\t")

write.table(merged[,c("targetid","Species")],file=paste0(paste0(projectdir,"/annotation/popmap_",projectname)),quote=F,row.names=F,col.names = F,sep="\t")
