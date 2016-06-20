#!/bin/bash

echo "SNP_Pipeline";
echo "====================";

#output parameters: should match those in FG/BG files_dir
compmodeOutput="BG" 
output="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/Output/"

#shared parameters:
sourcedir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/scripts/"
jobmode="SOCK"
nCPU=10

#createLDBlocks parameters:
inputsnplist="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/Background"
ldblockfile=$output$compmodeOutput"_createLDBlocks.block/Background"
ldblockdir=$output$compmodeOutput"_createLDBlocks.block/"
lddir="/homes/olymp/hoor.al-hasani/Projects/SNPS/Linkage_Data/LD_Hoor_liftOver/LD_hg19/LD_19/HapMap_files.RData"
loadhapmap="T"
joinBlocks="FALSE"
R=0.9

#randLDBlocks parameters:
Foreground="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/Output/FG_dupBlocks.block/Unique_FG_MSeQTLs.block"
plotTitle="Sampling_weights_Freq"
iter=100
plot="TRUE"
randblockdir=$output$compmodeOutput"_randLDBlocks.block/"

#dupBlocks
dupblockdir=$output$compmodeOutput"_dupBlocks.blockx/"
plotTitle_dup="Duplicated_blocks_vs._Rounds"
FG_mode="FALSE"
dupblockoutput=$dupblockdir


mkdir -p $ldblockdir
mkdir -p $randblockdir
mkdir -p $dupblockdir

echo "Output dir:"
echo $ldblockdir
echo $dupblockdir
echo $randblockdir
echo "====================";

#Background
#echo "createLDBlocks.R";
#time R --no-save --no-restore --slave --file=createLDBlocks.R --args $sourcedir $inputsnplist $lddir $ldblockdir $R $jobmode $nCPU $loadhapmap $joinBlocks 2> $output$compmodeOutput"_createLDBlocks.log" >$output$compmodeOutput"_createLDBlocks.output"
#echo "====================";

#Background sampling
#echo "randLDBlocks.R";
#time R --no-save --no-restore --slave --file=randLDBlocks.R --args $sourcedir $ldblockdir $randblockdir $iter $Foreground $jobmode $plot $plotTitle $nCPU 2> $output$compmodeOutput"_randLDBlocks.log" >$output$compmodeOutput"_randLDBlocks.output"
#echo "====================";

#Background filtering
echo "dupBlocks.R";
time R --no-save --no-restore --slave --file=dupBlocks.R --args $sourcedir $randblockdir $dupblockoutput $jobmode $nCPU $plotTitle_dup $FG_mode 2> $output$compmodeOutput"_dupBlocksx.log" >$output$compmodeOutput"_dupBlocksx.output"
echo "====================";

