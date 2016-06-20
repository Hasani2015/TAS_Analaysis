#!/bin/bash

echo "SNP_Pipeline";
echo "====================";

#output parameters: should match those in FG/BG files_dir
compmodeOutput="FG"
output="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/Output/New_filtering/TAS/Step_3/"

#shared parameters:
sourcedir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/scripts/"
jobmode="SOCK"
nCPU=10

#createLDBlocks parameters:
inputsnplist="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/Output/New_filtering/TAS/Step_2/"
lddir="/homes/olymp/hoor.al-hasani/Projects/SNPS/Linkage_Data/LD_Hoor_liftOver/LD_hg19/LD_19/HapMap_files.RData"
loadhapmap="T"
joinBlocks="F"
R=0.9
ldblockdir=$output$compmodeOutput"_createLDBlocks.block/"

# map the blocks - original TAS list
tas_dir="/homes/olymp/hoor.al-hasani/Projects/SNPS/output/TAS/TAS_Filter/"
tas_block_dir=$ldblockdir
output="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/Output/New_filtering/TAS/Step_4/"
map_dir=$output$compmodeOutput"_mapBlocks.map.block/"

#find TAS in platforms 
platforms_dir="/homes/olymp/hoor.al-hasani/Projects/SNPS/output/SNP_Output/Step_4/"
TAS_dir=$map_dir
output="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/Output/New_filtering/TAS/Step_5/"
tas_in_platforms_dir=$output$compmodeOutput"_TASinACS.map.block/"  
  
#weight tas
TAS_dir=$tas_in_platforms_dir
output="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/Output/New_filtering/TAS/Step_6/"
tas_w_out=$output$compmodeOutput"_TAS_W.block/TAS_W_new_filtering" 
tas_w_dir=$output$compmodeOutput"_TAS_W.block/"
joinBlocks="FALSE"
W_basis="Pubmed_ID"

#dupBlocks
#Foreground=$output$compmodeOutput"_createLDBlocks.block/TAS_fake_w.block"
plotTitle_dup="Duplicated_blocks_vs._Rounds"
FG_mode="TRUE"
dupblockdir=$output$compmodeOutput"_dupBlocks.block/"
dupblockoutput=$dupblockdir"Unique_FG.block"

mkdir -p $ldblockdir
mkdir -p $map_dir
mkdir -p $tas_in_platforms_dir
mkdir -p $tas_w_dir

#mkdir -p $dupblockdir


echo "Output dir:"
echo $ldblockdir
echo $map_dir
echo $tas_in_platforms_dir
echo $tas_w_dir
#echo $dupblockdir

#echo "====================";
#echo "createLDBlocks.R";
#time R --no-save --no-restore --slave --file=createLDBlocks.R --args $sourcedir $inputsnplist $lddir $Foreground $R $jobmode $nCPU $loadhapmap $joinBlocks 2> $output"createLDBlocks.log" >$output"createLDBlocks.output"
#echo "====================";

#echo "map_blocks_ACS.R";
#time R --no-save --no-restore --slave --file=map_blocks_ACS.R --args $tas_dir $tas_block_dir $map_dir $sourcedir 2> $output"mapBlocks.log" >$output"mapBlocks.output"

#echo "====================";
#echo "findTAS_Platforms.R";
#time R --no-save --no-restore --slave --file=findTAS_Platforms.R --args $sourcedir $platforms_dir $TAS_dir $tas_in_platforms_dir $jobmode $nCPU 2> $output"tas_acs.log" >$output"tas_acs.output" 
echo "====================";
echo "TAS_weight.R";
time R --no-save --no-restore --slave --file=TAS_weight.R --args $sourcedir $TAS_dir $tas_w_dir $jobmode $nCPU $joinBlocks $W_basis 2> $output"tas_w.log" >$output"tas_w.output"

#Foreground filtering
#echo "dupBlocks.R";
#time R --no-save --no-restore --slave --file=dupBlocks.R --args $sourcedir $Foreground $dupblockoutput $jobmode $nCPU $plotTitle_dup $FG_mode 2> $output"dupBlocks.log" >$output"dupBlocks.output"
#echo "====================";






