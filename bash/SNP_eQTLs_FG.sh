#!/bin/bash

echo "SNP_Pipeline";
echo "====================";

#output parameters: should match those in FG/BG files_dir
compmodeOutput="FG" #/FG/""
output="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/Output/"

#shared parameters:
sourcedir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/scripts/"
jobmode="SOCK"
nCPU=10

#createLDBlocks parameters:
inputsnplist="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/eqtl_snps"
ldblockfile=$output$compmodeOutput"_createLDBlocks.block/eqtl_snps"
ldblockdir=$output$compmodeOutput"_createLDBlocks.block/"
lddir="/homes/olymp/hoor.al-hasani/Projects/SNPS/Linkage_Data/LD_Hoor_liftOver/LD_hg19/LD_19/HapMap_files.RData"
loadhapmap="T"
joinBlocks="TRUE"
R=0.9

#dupBlocks
dupblockdir=$output$compmodeOutput"_dupBlocks.block/"
plotTitle_dup="Duplicated_blocks_vs._Rounds"
FG_mode="TRUE"
dupblockoutput=$dupblockdir"Unique_FG_MSeQTLs.block"


mkdir -p $ldblockdir
mkdir -p $dupblockdir

echo "Output dir:"
echo $ldblockdir
echo $dupblockdir
echo "====================";

#Foreground
echo "createLDBlocks.R";
time R --no-save --no-restore --slave --file=createLDBlocks.R --args $sourcedir $inputsnplist $lddir $ldblockfile $R $jobmode $nCPU $loadhapmap $joinBlocks 2> $output$compmodeOutput"_createLDBlocks.log" >$output$compmodeOutput"_createLDBlocks.output"
echo "====================";

#Foreground filtering
echo "dupBlocks.R";
time R --no-save --no-restore --slave --file=dupBlocks.R --args $sourcedir $ldblockfile $dupblockoutput $jobmode $nCPU $plotTitle_dup $FG_mode 2> $output$compmodeOutput"_dupBlocks.log" >$output$compmodeOutput"_dupBlocks.output"
echo "====================";

