######################################################################
# PURPOSE: create the LD_block
# INPUT: HapMap data frame of chromosome *, one row of SNP data frame 
# of the same chromosome
# RETURNS: data frame object corresponds to the given SNP
#####################################################################

snpp.ld_block <- function(df_Pos, LD){  
  
  ID = df_Pos["SNP_ID"]
  Pos = as.numeric(df_Pos["SNP_Pos"])
  W = as.numeric(df_Pos["SNP_W"])
  Chr = df_Pos["Chr"]
  
  snp = LD[LD$Pos1 == Pos | LD$Pos2 == Pos,]
  ldp_list = unique(c(snp$Pos1, snp$Pos2))
  
  if(length(ldp_list) >=1){
    sPos = min(ldp_list) 
    ePos = max(ldp_list)
  }else{
    sPos = Pos-1
    ePos = Pos
    ldp_list = Pos
  }
  
  print (paste("ID", ID))
  print (paste("Pos", Pos))
  print (paste("sPos", sPos))
  print (paste("ePos", ePos))
  print (paste("LDP", ldp_list))
  df_block = data.frame(SNP_ID = ID, SNP_Pos = Pos, SNP_W = W, Chr = Chr, sPos = sPos, ePos = ePos, LDP_list = paste(ldp_list, collapse = ","))

  return (df_block)
} # end function

##############################################################
# PURPOSE: filter HapMap data frame and call "snpp.ld_block"
#             one chromosome a time
# RETURNS: list of chromosome's blocks 
##############################################################

snpp.call_chr <- function(df_HapMap_chr, df_SNP_chr, R){
  
  colnames(df_HapMap_chr) = c("Pos1","Pos2","R")
  df_HapMap_chr = df_HapMap_chr[df_HapMap_chr$R > R, ]
  
  block_list = apply(df_SNP_chr,1, snpp.ld_block, LD = df_HapMap_chr)

  df_block = as.data.frame(matrix(unlist(block_list), ncol = 7))
  colnames(df_block) = c("SNP_ID","SNP_Pos", "SNP_W", "Chr", "sPos", "ePos", "LDP_list")
  print("=============")
  return (block_list)
}

##############################################################
# PURPOSE: sample over the blocks on SNP basis
#             one chromosome a time
# DEPENDS: plyr library
# RETURNS: list/file of chromosome's blocks of the given size 
##############################################################
snpp.sample_snps <- function(round, df_blocks, sample_size, output_dir){
  
  Round_i_sample = sample(df_blocks$SNP_Pos, sample_size, replace = F, prob = df_blocks$SNP_W)
  df_sample_snps = df_blocks[df_blocks$SNP_Pos %in% Round_i_sample,]  

  output_file_name = paste(output_dir,"Round_",round,"_Blocks_sampling.block", sep ="")
  snpp.saveOutput(df_sample_snps[1:sample_size,], output_file_name)  
  return (df_sample_snps)
}

##############################################################
# PURPOSE: sample over the blocks on SNP basis
#             one chromosome a time
# DEPENDS: plyr library
# RETURNS: list/file of chromosome's blocks of the given size 
##############################################################
count_snp_weight<- function(Round_df, weight_count_df){
  x = ddply(Round_df, .(SNP_Pos), "nrow")
  x$SNP_Pos = paste("Pos",x$SNP_Pos,sep="")
  weight_count_df[x$SNP_Pos] = x$nrow
  return  (weight_count_df)
}

################################################################
# PURPOSE: compute the mean of a subset in the given data.frame
#             according to the given parameter (cond)
# RETURNS: numeric value (mean) 
################################################################
snpp.get_mean<- function(cond, df, col){
  return  (mean(df[[col]][df[["Weights"]] == cond]))
}

#################################################################
# PURPOSE: call creatGI fo the data frames within the given list
# DEPENDS: snow library
# RETURNS: GI object list
################################################################
snpp.callDF <- function(CL, block_list, ldp, nCPU){
  if(length(block_list) < as.numeric(nCPU)){
	GI_list = lapply(block_list, snpp.createGI, ldp=ldp)
  }else{GI_list = parLapply(CL, block_list, snpp.createGI, ldp=ldp)}
  return(GI_list)
}

##############################################################
# PURPOSE: create GenomeInterval objects (GI) for a data.frame
# DEPENDS: genomeIntervals library
# RETURNS: GI object
##############################################################
snpp.createGI <- function(df_block, ldp){
  GI_object = GenomeIntervals(chromosome = df_block$Chr, start = df_block$sPos, end = df_block$ePos)
  if(ldp ==T){
    f = as.vector(df_block$LDP_list)
    GI_object@annotation$LDP = f
  }
  return(GI_object)
}

##############################################################
# PURPOSE: check the length of the data frames for BG_rounds
# RETURNS: True value only, otherwise the program stops.
##############################################################
snpp.check_myrounds <- function(outer_list){
  outer_list_element_length = lapply(outer_list, function(x){nrow(x)})
  print(outer_list_element_length)
  if(any(outer_list_element_length != outer_list_element_length[[1]])) {
    stop ("snpp.check_myrounds: one of the background list has unequal rounds size")
  }else{return (outer_list_element_length[[1]])}
}

##############################################################
# PURPOSE: check the total number of duplicated blocks in one
#                             round
# RETURNS: save the unique blocks
##############################################################
snpp.myDup_blocks <- function(block, output_dir){
  my_block <- subset(block, select = c("Chr","sPos","ePos"))
  my_savedBlock <- block[!(duplicated(my_block)),]
  snpp.saveOutput(my_savedBlock, output_dir)
  return(my_savedBlock)
  }

##############################################################
# PURPOSE: save one chromosome a time
# RETURNS: file for each chromosome
##############################################################
snpp.saveOutput <- function(df_chr, output_dir){
  write.table(df_chr, file = output_dir, quote = F, sep = "\t", col.names = T, row.names = F, eol ="\n")
}

#################################################################
# PURPOSE: take as an input a file with the name of the needed 
#     columns and the chr to add + adding weight column
# RETURNS: filtered table
#################################################################
snpp.Filter_df_Add_w_chr <- function(col_names, file_df, chr){
  col_names_vect = strsplit(col_names, ",")
  my_pos = grep(".*Pos.*", unlist(col_names_vect), perl =T, ignore.case = T)
  x = unlist(col_names_vect)
  x[my_pos]= "SNP_Pos"
  
  df = file_df[,unlist(col_names_vect)]
  df = transform(df, SNP_W = 1)
  df = transform(df, Chr = as.character(chr))
  colnames(df) = c(x,"SNP_W","Chr")
  
  return(df)
}

#################################################################
# PURPOSE: take merge two data.frames according to the chosen 
#       columns according one column in each
# RETURNS: new merged data.frame
#################################################################
snpp.merge_df <- function(df1, df2, bydf1, bydf2, colnames_df1, colnames_df2){
  colnames_df1_vect = strsplit(colnames_df1, ",")
  colnames_df2_vect = strsplit(colnames_df2, ",")
  my_df1_subset = df1[, unlist(colnames_df1_vect)]
  my_df2_subset = df2[, unlist(colnames_df2_vect)]
  
  df = merge(my_df1_subset, my_df2_subset, by.x = bydf1, by.y = bydf2, all = T)
  return(df)
}

##############################################################
# PURPOSE: find number of matchs in a data.frame to the 
#       given position
# RETURNS: data.frame
##############################################################
get_my_hitdf  <- function(tas_ldp, platforms_df){
  pos_hit_list = lapply(tas_ldp, function(pos){
    df_hit = platforms_df[platforms_df$Position == pos,]
    return(df_hit)
  })
}

##############################################################
# PURPOSE: find if the given tas exists in the platforms
# DEPENDS: plyr library, stringr
# RETURNS: the found tas or NULL
##############################################################
snpp.tas_in_platforms <- function(tas_row, platforms_df){

   tas_ldp = unlist(strsplit(as.character(tas_row["LDP_list"]), ","))
   tas_pos_vect = unique(tas_ldp)                                                     # + vectorize unique list of tas ldp [tas, is included so its blocks sPos,ePos]
   pos_hit_list = get_my_hitdf(tas_pos_vect, platforms_df)
   pos_hit_df = ldply(pos_hit_list, data.frame)
   l = nrow(pos_hit_df)
   if(is.integer(l)){
     snp_ACS_vect = unlist(strsplit(as.character(tas_row["ACS"]), ","))                 # + vectorize the associated platforms
     acs_list = c(as.character(pos_hit_df$Accession),as.character(pos_hit_df$Comment))
     unique_acs = unique(unlist(str_extract_all(string = acs_list, pattern="GPL\\d+")))
     
     valid_acs = snp_ACS_vect %in% unique_acs
     
     if(is.integer(length(snp_ACS_vect[valid_acs])) & length(snp_ACS_vect[valid_acs])>0) {
       #SNP_ID      Pos Pubmed_ID             ACS SNP_W   Chr     sPos     ePos                   LDP_list
       df = data.frame(SNP_ID = tas_row["SNP_ID"], SNP_Pos = tas_row["SNP_Pos"], Pubmed_ID = tas_row["Pubmed_ID"],  ACS = tas_row["ACS"], SNP_W = tas_row["SNP_W"], Chr = tas_row["Chr"], sPos = tas_row["sPos"], ePos = tas_row["ePos"], LDP_list = tas_row["LDP_list"])    
       return (df)}
     } 
}
##############################################################
# PURPOSE: parse the tas_df tas by tas (row by row)
# DEPENDS: plyr library
# RETURNS: list of found tas chromosome and -1 for not
##############################################################
snpp.call_tas_chr <- function(tas_chr, platforms_chr){

  found_tas_list = vector('list', nrow(tas_chr))
  found_tas_list = apply(tas_chr, 1, snpp.tas_in_platforms, platforms_df = platforms_chr) 
  
 return(found_tas_list)
}

