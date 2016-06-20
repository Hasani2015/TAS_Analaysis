# load libraries and user functions
#source("/data/bioinf/projects/code/projects_misc/Tools/R/UserFunctions/genomeIntervals_userFunctions.R")
#source("/data/bioinf/projects/code/projects_misc/Tools/R/UserFunctions/genomeIntervals_userFunctions_enrichments.R")  # wrapper functions for calculating and plotting the enrichments
rm(list = ls())
# get input
args <- commandArgs(trailingOnly=T)
if(length(args) > 0){
  source_dir      = args[1] #"/media/Hasani/SNP_Pipeline/"
  input_SNP_file  = args[2] #"/media/Hasani/SNP_Pipeline/SNP_List_e"
  HapMap_dir       = args[3] #"/media/Hasani/SNP_Pipeline/LD/ld_e/
  output_dir      = args[4] #"/media/Hasani/SNP_Pipeline/" if FG=F, else /media/Hasani/SNP_Pipeline/my_foreground_something
  r               = args[5] #0.9
  job.mode        = args[6]
  nCPU            = args[7]
  loadHapMap      = args[8] # T/F
  joinBlocks 	  = args[9] # T/F
} else{stop (paste("Missing input.", "argument number", length(args)))}

print("Your input:")
print("========================================")
print(paste("source_dir", source_dir))
print(paste("input_SNP_file", input_SNP_file)) 
print(paste("HapMap_dir", HapMap_dir))
print(paste("output_dir", output_dir))
print(paste("r", r))
print(paste("job.mode", job.mode))
print(paste("nCPU", nCPU))
print(paste("loadHapMap", loadHapMap))
print(paste("joinBlocks", joinBlocks))

#source_dir      = "/media/Hasani/SNP_Pipeline/"
#input_SNP_file  = "/media/Hasani/SNP_Pipeline/SNP_list.snp"
#HapMap_dir       = "/media/Hasani/SNP_Pipeline/LD/ld/"
#output_dir      = "/media/Hasani/SNP_Pipeline/"
#r               = 0.9
#job.mode        = "SOCK"
#loadHapMap      = F
#nCPU = 3
source(paste(source_dir,"Utlis.R",sep=""))
source(paste(source_dir,"SNPP_function.R",sep=""))

library(plyr)           # arrange and ddply needed for ggplot2 of ecdf
library(snow)           # for parallel version of sapply == parSapply

# initialize parallel version if needed
CL <- NULL
#nCPU=1
CL <- my.init.parallelRJob.snow(nCPU=nCPU, mode=job.mode, cl=CL)

clusterExport(CL, "source_dir")
clusterEvalQ(CL, library(plyr))   
clusterEvalQ(CL, source(paste(source_dir,"Utlis.R",sep="")))
clusterEvalQ(CL, source(paste(source_dir,"SNPP_function.R",sep="")))

# Read input file
# from: data/testInput.txt
print(system.time(list_snp_input <- snpp.readInput(input_SNP_file))) # one entry per chromosome
if(length(names(list_snp_input)) ==0 ) {
  stop(paste("createLDBlocks.R: SNP-chromosomes are missing!"))}

print("Done with the SNPs")

# Read HapMap data
# from: /data/bioinf/databases/proteom_databases/code/hapmap/v27/data/

if(loadHapMap){
	print("Loading LD_files")
	print(system.time(load(HapMap_dir)))
}else{
	print("Reading LD_files")
	print(system.time(
	list_hapmap <- snpp.readHapMap(CL, HapMap_dir))) # one entry per chromosome
      save(list_hapmap, file = paste(HapMap_dir,"/LD_files.RData", sep =""))
}
 
print("LD_files are done")
# all SNP-Chromosomes exist in HapMap
snp_chr_list = names(list_snp_input)
hm_chr_list = names(list_hapmap)
valid_chr = snp_chr_list %in% hm_chr_list

if(length(valid_chr) < length(names(list_snp_input)) ) {
        stop(paste("createLDBlocks.R: SNP-chromosomes don't exist in LD_files"))}

# The LD block calculation goes here
valid_HM_chr_list = names(list_hapmap) %in% names(list_snp_input)
print(names(list_hapmap))
print(names(list_snp_input))

print("=============================")
if(length(valid_HM_chr_list) != length(list_snp_input)) {
  stop(paste("createLDBlocks.R: SNP-chromosomes and LD_files mismatch!"))
}
vect_names = names(list_hapmap[valid_HM_chr_list])
list_blocks_list = clusterMap(CL, snpp.call_chr, list_hapmap[vect_names], df_SNP_chr = list_snp_input[vect_names], R =r)

df = parLapply(CL, list_blocks_list, ldply, data.frame)

output_list =  paste(output_dir, vect_names, ".block", sep ="")
if(joinBlocks==TRUE){
   joinedBlocks_df = ldply(df, data.frame)
   #output_file =  paste(output_dir, ".block", sep ="")
	#print(output_file)
   snpp.saveOutput(joinedBlocks_df, output_dir)
	print("DONE saving")

}else{clusterMap(CL, snpp.saveOutput, df, output_list)}

print("_______________________________________")

# save and close everything
#save.image(file="data/.RData")
stopCluster(CL)
