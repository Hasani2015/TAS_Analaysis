rm(list = ls())

args <- commandArgs(trailingOnly=T)
if(length(args) > 0){
  source_dir        = args[1] #"/media/Hasani/SNP_Pipeline/"
  platforms_dir     = args[2] #"/media/Hasani/SNP_Pipeline/BG_files_file"
  TAS_dir           = args[3] #"/media/Hasani/SNP_Pipeline/"
  output_dir        = args[4] #"/media/Hasani/SNP_Pipeline/"
  job.mode          = args[5]
  nCPU              = args[6]
} else{stop (paste("Missing input.", "argument number", length(args)))}

library(stringr)
library(plyr)          
library(snow)           # for parallel version of sapply == parSapply

options(warn = 1)
#source_dir        = "/media/Hasani/SNP_Pipeline/scripts/"
#platforms_dir     = "/media/Hasani/SNP_Pipeline/chnge_tas_filter.test/Step_4/"
#TAS_dir           = "/media/Hasani/SNP_Pipeline/chnge_tas_filter.test/"
#output_dir        = "/media/Hasani/SNP_Pipeline/chnge_tas_filter.test/filtered/"
#job.mode          = "SOCK"
#nCPU              = 3

print(paste("source_dir", source_dir))
print(paste("platforms_dir", platforms_dir)) 
print(paste("TAS_dir", TAS_dir))
print(paste("output_dir", output_dir))
print(paste("job.mode", job.mode))
print(paste("nCPU", nCPU)) 

source(paste(source_dir,"Utlis.R",sep=""))
source(paste(source_dir,"SNPP_function.R",sep=""))

# initialize parallel version if needed
CL <- NULL
CL <- my.init.parallelRJob.snow(nCPU=nCPU, mode=job.mode, cl=CL)

clusterExport(CL, "source_dir")
clusterEvalQ(CL, library(plyr))
clusterEvalQ(CL, library(stringr))   
clusterEvalQ(CL, source(paste(source_dir,"Utlis.R", sep="")))
clusterEvalQ(CL, source(paste(source_dir,"SNPP_function.R", sep="")))

print ("Reading TAS")
tas_list = snpp.readfiles.2Modes(CL = CL, dir_name = TAS_dir, reg_exp = ".*chr[\\dXY]*.*\\.map\\.block$",  N = 9 , nCPU = nCPU, Header = T, Fill=F)

# read platforms
print ("Reading Platforms (chromosome-wise)")
platforms_list = snpp.readfiles.2Modes(CL = CL, dir_name = platforms_dir, reg_exp = ".*chr.*[\\dXY]$",  N = 4 , nCPU = nCPU, Header = T, Fill = T)

if(length(tas_list) != length(platforms_list)){
  stop(paste("findTAS_Platforms.R: files mismatch! the given TAS and platforms files are not of the same length (missing chromosomes) "))}

print(paste("found tas files", length(tas_list)))
print(paste("found platforms files", length(platforms_list)))


# name both lists by chr
names(tas_list) =unlist(lapply(names(tas_list), function(my_names){
  l = gsub(".+chr","chr", my_names, perl =T,  ignore.case = T)
  l = paste("chr",gsub("\\D", "", l),sep ="")
  return (l)
}))

names(platforms_list)= unlist(lapply(names(platforms_list), function(my_names){
  l = gsub(".+chr","chr", my_names, perl =T,  ignore.case = T)
  l = paste("chr",gsub("\\D", "", l),sep ="")
  return (l)
}))

print("after renaming")

tas_names_vect = names(tas_list)
platforms_names_vect = names(platforms_list)

valid_chr = tas_names_vect[tas_names_vect %in% platforms_names_vect]

print("valid_chr")
print(valid_chr)

if(length(valid_chr) < length(tas_names_vect)){
  stop(paste("findTAS_Platforms.R: some TAS-chromosomes don't exist in platforms-chromosomes"))}

original_tas = lapply(tas_list[valid_chr], function(l) return(length(unique(l$SNP_Pos))))
names(original_tas) = valid_chr

# mapply tas_chr, platform_chr 
if(nCPU > length(tas_list[valid_chr])){
  x = mapply(snpp.call_tas_chr, tas_chr = tas_list[valid_chr], platforms_chr = platforms_list[valid_chr])
  df = lapply(x, ldply, data.frame)
  
}else{
  print("Paralleling over chromosomes")
  x = clusterMap(CL, snpp.call_tas_chr, tas_list[valid_chr], platforms_list[valid_chr])
  df = lapply(x, ldply, data.frame)
  }

# output snps_chr_i(1)
p = lapply(df, function(chr_df){
        noNA_df = chr_df[complete.cases(chr_df),]
        filtered_df = noNA_df[which(noNA_df$SNP_ID != -1),]
        return (filtered_df)
})

# save the output
final_tas = lapply(p, nrow)
names(final_tas) = valid_chr

print("Total number of original TAS: ")
print(paste(names(original_tas) , original_tas, sep = " = "))

print("Total number of filtered TAS: ")
print(paste(names(final_tas) , final_tas, sep = " = "))

my_df1 = t(do.call("rbind", original_tas))
my_df2 = t(do.call("rbind", final_tas))

my_df = rbind(my_df1,my_df2)
my_df = transform("file" = c("original TAS(unique)", "filtered TAS(not unique)"), my_df)

snpp.saveOutput (my_df, paste(output_dir,"/Read_ME_TAS_filtering_summary.csv", sep=""))
output_name_vect = paste(output_dir, "/TAS_",valid_chr, ".map.block", sep="")
save(p,file = paste(output_dir, "P.RData", sep="") )
f = mapply(snpp.saveOutput, p, output_name_vect)

# save and close everything
stopCluster(CL)
