

data1<-read.csv("samples/DGal21-6040_samples-Table\ 1_galaxids.csv")
data2<-read.csv("samples/targets_galaxiids_H2HL5DMXY_1.csv")
data3<-read.csv("samples/galaxiids_20210322_AWedited3Apr2021.csv")

merged<-merge(data1,data2,by.x="id",by.y="genotype",all.y=T)
merged<-merge(merged,data3,by.x="id",by.y="Genotype")

merged$Species[grepl("olidus",merged$pop)]<-merged$pop[grepl("olidus",merged$pop)]
merged$Species<-gsub("G_olidus_","Galaxias olidus ",merged$Species)



merged<-merged[,c("targetid","Species")]
write.table(merged,file="popmap_0702CR30_2",sep="\t",quote=F,row.names=F)