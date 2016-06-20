# load libraries and user functions
#source("/data/bioinf/projects/code/projects_misc/Tools/R/UserFunctions/genomeIntervals_userFunctions.R")
#source("/data/bioinf/projects/code/projects_misc/Tools/R/UserFunctions/genomeIntervals_userFunctions_enrichments.R")  # wrapper functions for calculating and plotting the enrichments
rm(list = ls())
# get input
args <- commandArgs(trailingOnly=T)
if(length(args) > 0){
  source_dir        = args[1]
  dir_background    = args[2] 
  output_dir        = args[3]
  job.mode          = args[4]
  nCPU              = args[5]
  plots_name        = args[6]
  FG_mode           = args[7]
} else{stop (paste("Missing input.", "argument number", length(args)))}

library(plyr)           # arrange and ddply needed for ldply and ddply
library(snow)           # for parallel version of sapply == parSapply
library(ggplot2)
options(warn = 1)

print(paste("source_dir", source_dir))
print(paste("dir_background", dir_background))
print(paste("output_dir", output_dir))
print(paste("job.mode", job.mode))
print(paste("nCPU", nCPU))
print(paste("plots_name", plots_name))
print(paste("FG_mode", FG_mode))

source(paste(source_dir,"Utlis.R",sep=""))
source(paste(source_dir,"SNPP_function.R",sep=""))
# initialize parallel version if needed
CL <- NULL
CL <- my.init.parallelRJob.snow(nCPU=nCPU, mode=job.mode, cl=CL)

clusterExport(CL, "source_dir")
clusterEvalQ(CL, library(plyr)) 
clusterEvalQ(CL, library(genomeIntervals))
clusterEvalQ(CL, source(paste(source_dir,"Utlis.R", sep="")))
clusterEvalQ(CL, source(paste(source_dir,"SNPP_function.R", sep="")))

# Read input file
# from: data/testInput.txt
if(FG_mode ==TRUE) {
  print("FG-Mode")
  df_fgBlocks = snpp.readFile(dir_background, T)
}else{
  print("BG-Mode")
  Background_list = lapply(dir_background, snpp.readBlocks, CL = CL, reg_exp = "^Round\\_\\d+.+\\.block$", N = 7, nCPU)
  check_myBG_roundsize = sapply(Background_list, snpp.check_myrounds, simplify = T)
  sample_size = check_myBG_roundsize
  print(sample_size)
}

# compute
if(FG_mode ==TRUE) {
  print("Warrning: in FG_mode the outputdir should also contain the name of the file to be saved (e.g. ~/my_new_unique_FG.block) ")
  my_blocks <- subset(df_fgBlocks, select = c("Chr","sPos","ePos"))
  my_dupblocks_nr = nrow(my_blocks[duplicated(my_blocks),])
  my_uniqBlocks = df_fgBlocks[!duplicated(my_blocks),]
  snpp.saveOutput(df_chr = my_uniqBlocks ,output_dir = output_dir)
}else{
  print("about to sapply my_dupBlocks_list")
 
  my_dupBlocks_list = data.frame(sapply(Background_list, function(my_df){
    f = parLapply(CL,my_df, function(block_df){
      my_block <- subset(block_df, select = c("Chr","sPos","ePos"))
      return (nrow(my_block[duplicated(my_block),]))})
    return(f)
  }))
 
my_dupBlocks_list = cbind(rownames(my_dupBlocks_list), my_dupBlocks_list)
 
  colnames(my_dupBlocks_list) = c("round","dupBlock")
  rownames(my_dupBlocks_list) = 1:nrow(my_dupBlocks_list)
 
  print(paste("Saving the unique blocks in", output_dir))
  f = lapply(Background_list, function(x){
    df = x[my_dupBlocks_list$round]
    names(df) = gsub(".+\\/Round","Round",names(df),perl =T)
    outputname = paste(output_dir, names(df), sep="")
    d = clusterMap(CL, snpp.myDup_blocks, block = df, output_dir = outputname)
    return(df)})
 
  print("DONE")
  print("Saving plot")
  library(ggplot2)
  my_dupBlocks_list$round = gsub(".+\\/Round","Round",my_dupBlocks_list$round,perl =T)
  nr.dubBlocks = as.numeric(my_dupBlocks_list$dupBlock)
  my_title = paste("Duplicated blocks vs. Rounds (sample_size=", sample_size,")")
  p =qplot(x = round, y = as.numeric(dupBlock), data = my_dupBlocks_list, geom = "point") + labs(title = my_title, x = "Round", y = "Total number of duplicated blocks of each rounds")+ theme(axis.text.x=element_text(angle=-90, hjust=0, size =5))
  p
  ggsave(paste(output_dir,plots_name,".pdf",sep=""), plot = p, width = 24)
}
# close cluster
stopCluster(CL)
