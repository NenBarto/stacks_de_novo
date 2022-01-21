#conda init 
#conda activate stacks

#this script runs stacks on the optimised parameter

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
for m in {3..3};do
		outDir=$stackDir"optimised"$m"_"$m
		mkdir -p $outDir

denovoLine="denovo_map.pl \
--samples $inDir \
--popmap $annotationDir/popmap_0702CR30 \
--out-path $outDir -T $ncores \
--min-samples-per-pop 0.9 \
-X \"populations: --vcf --genepop\" \
-M $m -n $m --resume" 
qsub -b y -hold_jid opt$m$n -wd $logDir -j y -N opt$m$m -R y -pe smp $ncores -V $denovoLine

done;








