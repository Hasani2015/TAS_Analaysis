rm(list = ls())
args <- commandArgs(trailingOnly=T)
if(length(args) > 0){
  file        = args[1] 
  col_names   = args[2] 
  chr         = args[3] 
  output_file = args[4]
  sourcedir   = args[5]
} else{stop (paste("Missing input.", "argument number", length(args)))}

print(paste("file", file))
print(paste("col_names", col_names))
print(paste("chr", chr))
print(paste("output_file", output_file))

#file = "/media/Hasani/SNPS/output/TAS/TAS_Filter/Chr_23_TAS_PubID_ACS"
#chr = "chr23"
#col_names = "SNP_ID,Pos"
#output_file = "/media/Hasani/SNPS/output/TAS/TAS_2Blocks_W1/Chr_23_TAS" 
#sourcedir = "/media/Hasani/SNP_Pipeline/scripts/"
#==============
source(paste(sourcedir,"SNPP_function.R", sep ="/"))
file_df = read.table(file, header = T)
my_df = snpp.Filter_df_Add_w_chr(col_names, file_df, chr)

write.table(my_df, output_file, quote=FALSE, row.names=FALSE, col.names=T, sep="\t")
