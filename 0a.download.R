library(rdrop2)

#token <- drop_auth(rdstoken="droptoken.rds")
#saveRDS(token, "droptoken.rds")

token <- readRDS("droptoken.rds")
# Then pass the token to each drop_ function
drop_acc(dtoken = token)

#inDir="/cesar Team Folder/_Customer Project Files/0702 DELWP (previously DSE)/0702CR30 Genomics of Mountain Galaxias/7. Data/raw data"
inDir="/test"

files<-drop_dir(inDir,dtoken = token)

#download all files
for(i in 1:dim(files)[1]){
	cat(i)
	drop_download(as.character(files[i,3]),dtoken = token)
}

inDir="/cesar Team Folder/_Customer Project Files/0702 DELWP (previously DSE)/0702CR27 Dingo scat project/7.) Data/Dartseq output-DRef20-5460"
files<-drop_dir(inDir,dtoken = token)

#download all files
for(i in 1:dim(files)[1]){
	cat(i)
	drop_download(as.character(files[i,3]),dtoken = token)
}

x <- drop_search(inDir,dtoken = token)
drop_download(dtoken = token, path = inDir)



#curl -X POST https://content.dropboxapi.com/2/files/download --header "Authorization: Bearer iTRjXf6DX4MAAAAAAAAAAdPDyWQBxoKYM-G_RS7hU-CrXQXgdq6xZ9IfuD91Jy3y" --header "Dropbox-API-Arg: {\"path\": \"/cesar Team Folder/_Customer Project Files/0702 DELWP (previously DSE)/0702CR27 Dingo scat project/7.) Data/Output-DCan21-6106"}""
#
outDir="/cesar Team Folder/_Customer Project Files/0702 DELWP (previously DSE)/0702CR27 Dingo scat project/7.) Data/Inhouse_stacks"
drop_upload("combined_annotation.txt",dtoken = token, path = outDir)
