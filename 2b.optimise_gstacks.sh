#conda init 
#conda activate gatk

#download:
#curl -X POST https://content.dropboxapi.com/2/files/download --header "Authorization: Bearer iTRjXf6DX4MAAAAAAAAAAdPDyWQBxoKYM-G_RS7hU-CrXQXgdq6xZ9IfuD91Jy3y" --header "Dropbox-API-Arg: {\"path\": \"/cesar Team Folder/_Customer Project Files/0702 DELWP (previously DSE)/0702CR27 Dingo scat project/7.) Data/Output-DCan21-6106"}""
#curl -X POST https://content.dropboxapi.com/2/files/download --header "Authorization: Bearer iTRjXf6DX4MAAAAAAAAAAdPDyWQBxoKYM-G_RS7hU-CrXQXgdq6xZ9IfuD91Jy3y" --header "Dropbox-API-Arg: {\"path\": \"test/2467376.FASTQ.gz\"}" -o "./annotation.tar.gz"

#https://www.dropbox.com/home/test/2467376.FASTQ.gz
#/_Customer Project Files/0702 DELWP (previously DSE)/0702CR27 Dingo scat project/7.) Data/Output-DCan21-6106

#to do 
#run one sample 2209339
#GATK
#R package

projectname="0702CR30"

inDir="/share/ScratchGeneral/nenbar/projects/canfam/results/$projectname.step1radtags"
resultsDir="/share/ScratchGeneral/nenbar/projects/canfam/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/canfam/annotation"
index=$annotationDir/"CanFam3.1"
logDir="/share/ScratchGeneral/nenbar/projects/canfam/scripts/logs"
mkdir -p $logDir

bwaDir="$resultsDir/"$projectname".bwa"
mkdir -p $bwaDir
stackDir="$resultsDir/"$projectname".opt.stacks"
mkdir -p $stackDir

ncores=15
for m in {1..3};do
		outDir=$stackDir$m"_"$m
		mkdir -p $outDir

denovoLine="denovo_map.pl \
--samples $inDir \
--popmap $annotationDir/popmapOptimise_0702CR30 \
--out-path $outDir -T $ncores \
--min-samples-per-pop 0.8 \
-M $m -n $m --resume" 
qsub -b y -hold_jid opt$m$n -wd $logDir -j y -N opt$m$m -R y -pe smp $ncores -V $denovoLine

done;








