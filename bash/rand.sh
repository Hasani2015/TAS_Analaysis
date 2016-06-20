#!/bin/bash

echo "SNP_Pipeline";
echo "====================";

#output parameters:
compmodeOutput="" #"BG_" #/FG_/""
output="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/test_dup/"
ldblockdir=$output$compmodeOutput"createLDBlocks.block/"
randblockdir=$output$compmodeOutput"randLDBlocks.block/"
ovGIdir=$output$compmodeOutput"overlapGI1/"

#shared parameters:
sourcedir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/scripts/"
jobmode="SOCK"
nCPU=1

#randLDBlocks parameters:
Foreground="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/Foreground_list.df.block"
plotTitle="Sampling_weights_Freq"
iter=100
plot="F"
# + ldblockdir

#mkdir -p $ldblockdir
mkdir -p $randblockdir
#mkdir -p $ovGIdir

echo "Output dir:"
echo $ldblockdir
echo $randblockdir
echo $ovGIdir
echo "====================";

echo "randLDBlocks.R";
time R --no-save --no-restore --slave --file=randLDBlocks.R --args $sourcedir $ldblockdir $randblockdir $iter $Foreground $jobmode $plot $plotTitle $nCPU 2> $sourcedir"randLDBlocks.log" >$sourcedir"randLDBlocks.output"
echo "====================";
