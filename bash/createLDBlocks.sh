#!/bin/bash

echo "SNP_Pipeline";
echo "====================";

#output parameters: should match those in FG/BG files_dir
compmodeOutput="Imp"
output="/homes/olymp/hoor.al-hasani/Projects/SNPS/output/TAS/TAS_filtering_steps_source_files/Output/"

#shared parameters:
sourcedir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/scripts/"
jobmode="SOCK"
nCPU=5

#createLDBlocks parameters:
inputsnplist="/homes/olymp/hoor.al-hasani/Projects/SNPS/output/TAS/TAS_filtering_steps_source_files/TAS_byImp4LD"
lddir="/homes/olymp/hoor.al-hasani/Projects/SNPS/Linkage_Data/LD_Hoor_liftOver/LD_hg19/LD_19/HapMap_files.RData"
loadhapmap="T"
joinBlocks="T"
R=0.9
ldblockdir=$output$compmodeOutput"_createLDBlocks.block/"

mkdir -p $ldblockdir

echo "Output dir:"
echo $ldblockdir

echo "====================";
echo "createLDBlocks.R";
time R --no-save --no-restore --slave --file=createLDBlocks.R --args $sourcedir $inputsnplist $lddir $ldblockdir $R $jobmode $nCPU $loadhapmap $joinBlocks 2> $output"createLDBlocksimp.log" >$output"createLDBlocksimp.output"
echo "====================";




