#!/bin/bash

echo "SNP_Pipeline";
echo "====================";

#output parameters: should match those in FG/BG files_dir
compmodeOutput="FG" #/FG/""
output="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/Output/"
#ldblockdir=$output$compmodeOutput"_createLDBlocks.block/"
#randblockdir=$output$compmodeOutput"_randLDBlocks.block/"
dupblockdir=$output$compmodeOutput"_dupBlocks.block/"
#ovGIdir=$output$compmodeOutput"_overlapGI/"

#shared parameters:
sourcedir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/scripts/"
jobmode="SOCK"
nCPU=4

#createLDBlocks parameters:
inputsnplist="/homes/olymp/hoor.al-hasani/Projects/M.S.eQTLs/eqtl_snps"
#lddir="/homes/olymp/hoor.al-hasani/Projects/SNPS/Linkage_Data/LD_Holger/3col/"

#inputsnplist="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/SNP_List"
lddir="/homes/olymp/hoor.al-hasani/Projects/SNPS/Linkage_Data/LD_Hoor_liftOver/LD_hg19/LD_19/HapMap_files.RData"

loadhapmap="T"
joinBlocks="T"
R=0.9
#fg
ldblockfile=$output$compmodeOutput"_createLDBlocks.block/eqtl_snps"
#bg
#ldblockdir=$output$compmodeOutput"_createLDBlocks.block/"

#randLDBlocks parameters:
#Foreground=$inputsnplist".block"
#plotTitle="Sampling_weights_Freq"
#iter=100
#plot="TRUE"
# + ldblockdir

#dupBlocks
plotTitle_dup="Duplicated_blocks_vs._Rounds"
FG_mode="TRUE"
dupblockoutput=$dupblockdir"Unique_FG_MSeQTLs.block"

#overlapGI parameters:
#BG_files_dir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/BG_files_file"
#FG_files_dir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/FG_files_file"
#AS_files_dir="/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/AS_files_file"
#compmode=2

mkdir -p $ldblockdir
#mkdir -p $randblockdir
mkdir -p $dupblockdir
#mkdir -p $ovGIdir

echo "Output dir:"
echo $ldblockdir
#echo $randblockdir
echo $dupblockdir
#echo $ovGIdir
echo "====================";

#Foreground
echo "createLDBlocks.R";
time R --no-save --no-restore --slave --file=createLDBlocks.R --args $sourcedir $inputsnplist $lddir $ldblockfile $R $jobmode $nCPU $loadhapmap $joinBlocks 2> $output$compmodeOutput"_createLDBlocks.log" >$output$compmodeOutput"_createLDBlocks.output"
echo "====================";

#Foreground filtering
#echo "dupBlocks.R";
#time R --no-save --no-restore --slave --file=dupBlocks.R --args $sourcedir $ldblockdir $dupblockoutput $jobmode $nCPU $plotTitle_dup $FG_mode 2> $output$compmodeOutput"_dupBlocks.log" >$output$compmodeOutput"_dupBlocks.output"
#echo "====================";

#echo "randLDBlocks.R";
#time R --no-save --no-restore --slave --file=randLDBlocks.R --args $sourcedir $ldblockdir $randblockdir $iter $dupblockoutput $jobmode $plot $plotTitle $nCPU 2> $sourcedir"randLDBlocks.log" >$sourcedir"randLDBlocks.output"
#echo "====================";

#Background filtering
#echo "dupBlocks.R";
#FG_mode="F"
 
#time R --no-save --no-restore --slave --file=dupBlocks.R --args $sourcedir $randblockdir $dupblockdir $jobmode $nCPU $plotTitle_dup $FG_mode 2> $sourcedir"dupBlocks.log" >$sourcedir"dupBlocks.output"
#echo "====================";


#ovGIdirAS=$ovGIdir"Exon/"
#mkdir -p $ovGIdirAS

#echo "overlapGI.R";
#echo "Exon"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIEx.log" >$sourcedir"overlapGIEx.output"
#echo "====================";
#ovGIdirAS=$ovGIdir"Incase_you_forgot_2create_output_dir/"
#mkdir -p $ovGIdirAS
#ovGIdirAS=$ovGIdir"Intron/"
#mkdir -p $ovGIdirAS

#echo "Intron"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIIn.log" >$sourcedir"overlapGIIn.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"3UTR/"
#mkdir -p $ovGIdirAS

#echo "3'UTR"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGI3U.log" >$sourcedir"overlapGI3U.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"5UTR/"
#mkdir -p $ovGIdirAS

#echo "5'UTR"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGI5U.log" >$sourcedir"overlapGI5U.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"Gene/"
#mkdir -p $ovGIdirAS

#echo "Gene"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIGe.log" >$sourcedir"overlapGIGe.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"MCS/"
#mkdir -p $ovGIdirAS

#echo "MCS"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIMCS.log" >$sourcedir"overlapGIMCS.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"CpG/"
#mkdir -p $ovGIdirAS

#echo "CpG_Islands"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGICpG.log" >$sourcedir"overlapGICpG.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"CC_intergenic/"
#mkdir -p $ovGIdirAS

#echo "CC_intergenic"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGICCintg.log" >$sourcedir"overlapGICCintg.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"CC_intronic/"
#mkdir -p $ovGIdirAS

#echo "CC_intronic"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGICCint.log" >$sourcedir"overlapGICCint.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"ECS_RNAz/"
#mkdir -p $ovGIdirAS

#echo "ECS_RNAz"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGICCECS_RNAz.log" >$sourcedir"overlapGIECS_RNAz.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"ECS/"
#mkdir -p $ovGIdirAS

#echo "ECS"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGICCECS.log" >$sourcedir"overlapGIECS.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"cars/"
#mkdir -p $ovGIdirAS

#echo "cars"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIcars.log" >$sourcedir"overlapGIcars.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"p53/"
#mkdir -p $ovGIdirAS

#echo "p53"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIp53.log" >$sourcedir"overlapGIp53.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"Stat3/"
#mkdir -p $ovGIdirAS

#echo "Stat3"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIStat3.log" >$sourcedir"overlapGIStat3.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"lnc+TUCP/"
#mkdir -p $ovGIdirAS

#echo "lnc+TUCP"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIlnc+TUCP.log" >$sourcedir"overlapGIlnc+TUCP.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"lncRNA/"
#mkdir -p $ovGIdirAS

#echo "lncRNA"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIlncRNA.log" >$sourcedir"overlapGIlncRNA.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"lincRNA/"
#mkdir -p $ovGIdirAS

#echo "lincRNA"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIlincRNA.log" >$sourcedir"overlapGIlincRNA.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"TUCP/"
#mkdir -p $ovGIdirAS

#echo "TUCP"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGITUCP.log" >$sourcedir"overlapGITUCP.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"gencode_trans_fltrd/"
#mkdir -p $ovGIdirAS

#echo "gencode_trans_fltrd"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIgncd_trans.log" >$sourcedir"overlapGIgncd_trans.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"gencode_exon_fltrd/"
#mkdir -p $ovGIdirAS

#echo "gencode_exon_fltrd"
##time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIgncd_exon.log" >$sourcedir"overlapGIgncd_exon.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"gencode_psdgnG/"
#mkdir -p $ovGIdirAS

#echo "gencode_psdgnG"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIpsdgnG.log" >$sourcedir"overlapGIpsdgnG.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"OregAnno/"
#mkdir -p $ovGIdirAS

#echo "OregAnno"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIOregAnno.log" >$sourcedir"overlapGIOregAnno.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"gencode_CDS/"
#mkdir -p $ovGIdirAS

#echo "gencode_CDS"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIgencode_CDS.log" >$sourcedir"overlapGIgencode_CDS.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"5UTR/"
#mkdir -p $ovGIdirAS

#echo "5UTR"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGI5UTR.log" >$sourcedir"overlapGI5UTR.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"3UTR/"
#mkdir -p $ovGIdirAS

#echo "3UTR"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGI3UTR.log" >$sourcedir"overlapGI3UTR.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"gencode.intergenic/"
#mkdir -p $ovGIdirAS

#echo "gncd_intrgnc"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIgncd_intrgnc.log" >$sourcedir"overlapGIgncd_intrgnc.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"sno-miRNA_Gene/"
#mkdir -p $ovGIdirAS

#echo "sno-miRNA_Gene"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIsno-miRNA_Gene.log" >$sourcedir"overlapGIsno-miRNA_Gene.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"sno-miRNA_3UTR/"
#mkdir -p $ovGIdirAS

#echo "sno-miRNA_3UTR"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIsno-miRNA_3UTR.log" >$sourcedir"overlapGIsno-miRNA_3UTR.output"
#echo "====================";

#ovGIdirAS=$ovGIdir"sno-miRNA_5UTR/"
#mkdir -p $ovGIdirAS

#echo "sno-miRNA_5UTR"
#time R --no-save --no-restore --slave --file=overlapGI.R --args $sourcedir $BG_files_dir $FG_files_dir $AS_files_dir $ovGIdirAS $jobmode $nCPU 2> $sourcedir"overlapGIsno-miRNA_5UTR.log" >$sourcedir"overlapGIsno-miRNA_5UTR.output"
#echo "====================";


