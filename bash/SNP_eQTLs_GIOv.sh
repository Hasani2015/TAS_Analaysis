#!/bin/bash

echo "SNP_Pipeline";
echo "====================";

#output parameters: should match those in FG/BG files_dir
compmodeOutput="GI" 
output="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/Output/"

#shared parameters:
sourcedir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/scripts/"
jobmode="SOCK"
nCPU=10

#overlapGI parameters:
BG_files_dir="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/BG_files_file"
FG_files_dir="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/FG_files_file.txt"
AS_files_dir="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/AS_files_file"
ovGIdir=$output$compmodeOutput"_overlapGI/"

mkdir -p $ovGIdir

echo "Output dir:"
echo $ovGIdir
echo "====================";
echo "overlapGI.R"

time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdir $jobmode $nCPU 2> $output$compmodeOutput"_overlapGI.log" >$output$compmodeOutput"_overlapGI.output"
echo "====================";

