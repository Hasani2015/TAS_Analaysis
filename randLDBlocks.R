# load libraries and user functions
#source("/data/bioinf/projects/code/projects_misc/Tools/R/UserFunctions/genomeIntervals_userFunctions.R")
#source("/data/bioinf/projects/code/projects_misc/Tools/R/UserFunctions/genomeIntervals_userFunctions_enrichments.R")  # wrapper functions for calculating and plotting the enrichments
rm(list = ls())
# get input
args <- commandArgs(trailingOnly=T)
if(length(args) > 0){
  source_dir        = args[1] #"/media/Hasani/SNP_Pipeline/"
  input_dir_blocks  = args[2] #"/media/Hasani/SNP_Pipeline/"
  output_dir        = args[3] #"/media/Hasani/SNP_Pipeline/"
  iter              = args[4] 
  FG_file           = args[5] 
  job.mode          = args[6]
  plots             = args[7]
  plots_name        = args[8]
  nCPU       	    = args[9]
} else{stop (paste("Missing input.", "argument number", length(args)))}

library(plyr)           # arrange and ddply needed for ldply and ddply
library(snow)           # for parallel version of sapply == parSapply

print(paste("source_dir", source_dir))
print(paste("input_dir_blocks", input_dir_blocks))
print(paste("output_dir", output_dir))
print(paste("iter", iter))
print(paste("FG_file", FG_file))
print(paste("job.mode", job.mode))
print(paste("plots", plots))
print(paste("plots_name", plots_name))
print(paste("nCPU", nCPU))

source(paste(source_dir,"Utlis.R",sep=""))
source(paste(source_dir,"SNPP_function.R",sep=""))

# initialize parallel version if needed
CL <- NULL
CL <- my.init.parallelRJob.snow(nCPU=nCPU, mode=job.mode, cl=CL)

clusterExport(CL, "source_dir")
clusterEvalQ(CL, library(plyr))   
clusterEvalQ(CL, source(paste(source_dir,"Utlis.R",sep="")))
clusterEvalQ(CL, source(paste(source_dir,"SNPP_function.R",sep="")))

# Read input file
# from: data/testInput.txt
#print("Reading the blocks")

blocks_list =  snpp.readBlocks(CL,input_dir_blocks, "^chr[\\dXY]*\\.block$", 8, nCPU)
blocks_df = ldply(blocks_list, data.frame)
drop = ".id"
blocks_df = blocks_df[!(names(blocks_df) %in% drop)]

#read foreground list for sample_size information
#print("Reading FG")
FG_df = try_read_table_function(FG_file, H=T)
if(class(FG_df) =="data.frame"){
  sample_size =nrow(FG_df)
}else{
  stop("randLDBlocks.R: invalid foreground list! should has a block-format")
}
print(paste("sample_size", sample_size))
# The sampling goes here

print("Start sampling...")
rounds_vector = c(1:iter)
rand_blocks_list = vector("list", length(blocks_list))#rounds_vector
system.time(rand_blocks_list <- parLapply(CL,rounds_vector, snpp.sample_snps, df_blocks = blocks_df, sample_size = sample_size, output_dir = output_dir))
print("_______________________________________")

# The plot goes here
if(plots==TRUE){
  library(ggplot2)
  rounds_vector = paste("Round", rounds_vector, sep ="")
  names(rand_blocks_list) = rounds_vector

  df = ldply(rand_blocks_list,data.frame)
# df = do.call("rbind", rand_blocks_list)
  SNPs_weight_range = unique(blocks_df$SNP_W)
  
  SNPs_Pos_vect = unique(blocks_df$SNP_Pos)  
  snp_pos_df = unique(data.frame(blocks_df[,c("SNP_Pos","SNP_W")]))
  df_snps_p = data.frame(rbind(replicate(length(SNPs_Pos_vect), 0)))
  colnames(df_snps_p) = paste("Pos",SNPs_Pos_vect,sep="")
  
  df_snps_p = as.data.frame(sapply(snp_pos_df$SNP_Pos, function(pos){
		return(length(df$SNP_Pos[df$SNP_Pos== as.numeric(pos)])/as.numeric(iter))
		}))
  snp_pos_df = cbind(snp_pos_df, df_snps_p)
  colnames(snp_pos_df) = c("Pos","W","Freq")
  
  SNPs_weight_r = as.data.frame(sapply(SNPs_weight_range, function(x){return(sum(snp_pos_df$Freq[snp_pos_df$W==x])/length(snp_pos_df$Freq[snp_pos_df$W==x]))}))

  SNPs_weight_r = cbind(SNPs_weight_range, SNPs_weight_r)
  colnames(SNPs_weight_r) = c("SNPs_weights","Freq")

  g = ggplot()  + geom_point(aes(SNPs_weights, Freq), data = SNPs_weight_r, fill = "steelblue", col = "darkblue") + scale_x_continuous(breaks=SNPs_weight_range)#+ theme(axis.text.x  = element_text(angle=90, hjust=1.2, size=9)) 

  output_file_name = paste(output_dir,paste(plots_name,"_",as.numeric(iter),".csv", sep=""), sep ="")
  snpp.saveOutput(SNPs_weight_r,paste(output_dir,plots_name, sep =""))
  output_file_name = paste(output_dir,paste(plots_name,"_",as.numeric(iter),".pdf", sep=""), sep ="")
  ggsave(output_file_name, plot = g, width = 24)
  
  }
# stop cluster
stopCluster(CL)
