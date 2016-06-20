for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22;
do
	time R --no-save --no-restore --slave --file=filter_file.R --args ~/Projects/SNPS/output/TAS/TAS_Filter/Chr_${i}_TAS_PubID_ACS "SNP_ID,Pos" "chr"${i} ~/Projects/SNP_Pipeline/Output/New_filtering/TAS/Step_2/Chr${i}_TAS_conv ~/Projects/SNP_Pipeline/scripts/
done;
