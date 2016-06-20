rm(list = ls())
args <- commandArgs(trailingOnly=T)
if(length(args) > 0){
  tas_file    = args[1] 
  tas_block   = args[2] 
  output_dir = args[3]
  sourcedir   = args[4]
} else{stop (paste("Missing input.", "argument number", length(args)))}

# change the script to read more than file

print(paste("tas_file", tas_file))
print(paste("tas_block", tas_block))
print(paste("output_dir", output_dir))
print(paste("sourcedir", sourcedir))

#tas_file = "/media/Hasani/SNPS/output/TAS/TAS_Filter_incon_gwas/" 
#tas_block = "/media/Hasani/SNP_Pipeline/Output/Output/FG/createLDBlocks.block/" 
#sourcedir = "/media/Hasani/SNP_Pipeline/scripts/"
#output_dir = "/media/Hasani/SNP_Pipeline/chnge_tas_filter.test/"

library(plyr)
source(paste(sourcedir,"SNPP_function.R", sep ="/"))
source(paste(sourcedir,"Utlis.R", sep ="/"))

tas_list = snpp.readfiles.2Modes(NULL, tas_file, "Chr\\_*[\\dXY]*\\_*TAS\\_PubID\\_ACS", 4, 0, Header = T, Fill=F)
tas_block_list =  snpp.readfiles.2Modes(NULL,tas_block, "^chr[\\dXY]*\\.block$", 8, 0, Header =T, Fill=F)

names(tas_list) =unlist(lapply(names(tas_list), function(my_names){
  l = gsub(".+chr","chr", my_names, perl =T,  ignore.case = T)
  l = paste("chr",gsub("\\D", "", l),sep ="")
  return (l)
}))

names(tas_block_list) =unlist(lapply(names(tas_block_list), function(my_names){
  l = gsub(".+chr","chr", my_names, perl =T,  ignore.case = T)
  l = paste("chr",gsub("\\D", "", l),sep ="")
  return (l)
}))

drop = ".id"
tas_block_list = lapply(tas_block_list, function(l){
  return(l[!colnames(l) %in% drop])
})
 
tas_names_vect = names(tas_list)
tas_block_names_vect = names(tas_block_list)

valid_chr = tas_names_vect[tas_names_vect %in% tas_block_names_vect]

if(length(valid_chr) < length(tas_names_vect)){
  print(paste("Warrning: map_blocks_ACS.R: lists length mismatch!"))}

#print("................................")
#print("Number of unique entery in blocks")
chr_nr_table = lapply(tas_block_list[valid_chr], function(l){return(length(unique(l$SNP_Pos)))})
chr_nr_table = do.call("rbind", chr_nr_table)
chr_nr_table = t(chr_nr_table)
x = snpp.saveOutput(chr_nr_table, paste(output_dir,"Blocks_original_unique.csv", sep =""))

#print("Number of unique entery in TAS list")
chr_nr_table = sapply(tas_list[valid_chr], function(l){return(length(unique(l$Pos)))})
chr_nr_table = t(chr_nr_table)
x = snpp.saveOutput(chr_nr_table, paste(output_dir,"TAS_original_unique.csv", sep =""))

# find tas in blocks
tas_pos_in_blocks_list = vector("list", length(tas_list[valid_chr]))

tas_pos_in_blocks_list = mapply(function(df1,df2){
  positions_vect = unique(df2$SNP_Pos)
  hit_list = lapply(positions_vect, function(pos){
    return(df1[df1$Pos == pos,])
  })
  return(hit_list)
}, tas_list[valid_chr], tas_block_list[valid_chr])

tas_pos_in_blocks_df_list = lapply(tas_pos_in_blocks_list, ldply, data.frame)

#print("................................")
#print("Number of unique TAS in blocks")
chr_nr_table = sapply(tas_pos_in_blocks_df_list, function(l){return(length(unique(l$Pos)))})
chr_nr_table = t(chr_nr_table)
x = snpp.saveOutput(chr_nr_table, paste(output_dir,"TAS_found_in_blocks_unique.csv", sep =""))


tas_blocks_map_list = mapply(function(dfb,dft){
  positions_vect = unique(dfb$SNP_Pos)
  
  mapped_df = lapply(positions_vect, function(pos){
    block_hit_df = dfb[dfb$SNP_Pos==pos,]
    block_hit_row = block_hit_df[1,]
    tas_hit_df = dft[dft$Pos==pos,]
    
    df = data.frame(SNP_ID = rep(block_hit_row$SNP_ID, nrow(tas_hit_df)), SNP_Pos = rep(block_hit_row$SNP_Pos, nrow(tas_hit_df)), Pubmed_ID = tas_hit_df$Pubmed_ID, ACS = tas_hit_df$ACS, SNP_W= rep(block_hit_row$SNP_W, nrow(tas_hit_df)), Chr= rep(block_hit_row$Chr, nrow(tas_hit_df)), sPos= rep(block_hit_row$sPos, nrow(tas_hit_df)), ePos = rep(block_hit_row$ePos, nrow(tas_hit_df)), LDP_list = rep(block_hit_row$LDP_list, nrow(tas_hit_df)))
    return (df)
  })
 return(mapped_df) 
}, tas_block_list[valid_chr], tas_pos_in_blocks_df_list)

tas_blocks_map_df_list = lapply(tas_blocks_map_list, ldply, data.frame)

#print("................................")
#print("Number of unique maped blocks-TAS list")
chr_nr_table = sapply(tas_blocks_map_df_list, function(l){return(length(unique(l$SNP_Pos)))})
chr_nr_table = t(chr_nr_table)
x = snpp.saveOutput(chr_nr_table, paste(output_dir,"blocks-TAS_unique.csv", sep =""))

output_file = paste(output_dir,valid_chr,".map.block", sep ="")
x = mapply(snpp.saveOutput, tas_blocks_map_df_list, output_file)
print("_______________________________________")
