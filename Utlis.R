######################################################
# PURPOSE: Check if clusters for parallel processing are already created by makeCluster
#          If mode==EVE:  send job to EVE cluster --> MPI takes care of parallel process --> cluster object "cl" already exists
#          If mode==SOCK: job is send to CPUs on local machine (== multi-threading) --> create cluster object "cl"
# DEPENDS: snow library
# RETURNS: An object "cl" containing the clusters
######################################################
my.init.parallelRJob.snow <- function(nCPU=1, mode="SOCK", cl=NULL){
  
  if(mode=="SOCK"){
    cl <- makeCluster(rep("localhost", nCPU), type = "SOCK")
  }
  else if(!is.null(cl) && mode=="EVE"){
    cl <- cl
  }
  else if(is.null(cl) && mode=="EVE"){
    # If MPI cluster is NULL try to reset it:
    cl <- getMPIcluster()
    if(is.null(cl)){
      stop("No MPI cluster found and reinitialization of MPI cluster failed!")
    }
  }
  else{
    stop(paste("No valid mode for parallel processing: ", mode, sep=""))
  }
  
  if(is.null(cl)){
    stop("Could not create a cluster object!")
  }
  
  return(cl)
}

######################################################
# PURPOSE: read a table format file
#          If H == T, read with header, else read without header
# RETURNS: either a file or -1
######################################################
try_read_table_function <- function(file, H){

f = try(read.table(file, header = H))  
if(class(f) != "try-error"){
  return (f)
}else{return (-1)}
}

###################################################################
# PURPOSE: check the format of the files
# RETURNS: through an error if wrong format was found
##################################################################
check_files <- function (list_chr, chr, errr_msg, nr)
{
  if(length(list_chr[[chr]]) != nr){
    stop(paste("check_files( ",errr_msg," ): invalid number of columns (",length(list_chr[[chr]]),") in ", chr, sep=""))
  }  
}

###################################################################
# PURPOSE: read HapMap files and create list of their data frames
#          If indx
# DEPENDS: snow library
# RETURNS: list of data frames (chromosomes) named by "chr*"
##################################################################
snpp.readHapMap <- function(CL, dir_name){
  
files_dir_list = dir(dir_name, full.names =T)
indx = grep("LD\\d*\\_chr[\\dXYM]*$", files_dir_list, perl = T,  ignore.case = T)
LD_files_list = parLapply(CL, files_dir_list[indx], try_read_table_function, H = F)

names(LD_files_list) = files_dir_list[indx]
names(LD_files_list) = gsub(".+\\_chr","chr", names(LD_files_list), perl =T)
parLapply(CL, names(LD_files_list), check_files, list_chr = LD_files_list, errr_msg ="HapMap", nr =3)

return (LD_files_list);
}

################################################################
# PURPOSE: read a table format file corresponding to the blocks
# RETURNS: list of data frames corresponding to the chromosomes
################################################################
snpp.readBlocks <- function(CL,dir_name, reg_exp, N, nCPU){
  
  files_dir_list = dir(dir_name, full.names = FALSE)  

  #indx = grep("^chr[\\dXY]*\\.block$", files_dir_list, perl = T,  ignore.case = T)
  indx = grep(reg_exp, files_dir_list, perl = T,  ignore.case = T)

  files_dir_list[indx] = paste(dir_name, files_dir_list[indx], sep ="")  
  if(length(files_dir_list[indx]) < as.numeric(nCPU)){
	LD_blocks_list = lapply(files_dir_list[indx], try_read_table_function, H =T)
  }else{LD_blocks_list = parLapply(CL,files_dir_list[indx], try_read_table_function, H =T)}
  
  names(LD_blocks_list) = files_dir_list[indx]
  names(LD_blocks_list) = gsub(".+chr","chr", names(LD_blocks_list), perl =T)
  parLapply(CL, names(LD_blocks_list), check_files, list_chr = LD_blocks_list, errr_msg ="SNPs_LD_Blocks", nr =N)
  return (LD_blocks_list);
}

##############################################################
# PURPOSE: read one single file
# RETURNS:  data frames (one file) 
##############################################################
snpp.readFile <- function(file_name, H){
  
  Annot_file = try(read.table(file_name, header = H))
  if(class( Annot_file) != "try-error"){
    return ( Annot_file)
  }else{stop(paste("snpp.readFile: file [", Annot_file, "] has invalid format", sep=""))}
}

######################################################
# PURPOSE: read a table format file
#          If H == T, read with header, else read without header
# DEPENDS: snow library
# RETURNS: An object "cl" containing the clusters
######################################################
snpp.readInput <- function(file_name){
  
file_fields_number = max(count.fields(file_name))
if(file_fields_number != 4){
        stop(paste("snpp.readInput: program terminated: invalid columns number"))}

df_SNP = try(read.table(file_name, header = T))
if(class(df_SNP) == "try-error"){
        stop(paste("snpp.readInput: program terminated: invalid format"))}

default_snps_colnames = c("SNP_ID","SNP_Pos","SNP_W","Chr")
if(any(!(default_snps_colnames %in% colnames(df_SNP)))){
  stop(paste("snpp.readInput: invalid colnames; default are SNP_ID,SNP_Pos,SNP_W, and Chr!"))
}

list_snp_chr = as.list(levels(df_SNP$Chr))
names(list_snp_chr) = gsub("chromosome","chr",levels(df_SNP$Chr), perl =T)

for (i in levels(df_SNP$Chr)){
        list_snp_chr[[i]] = df_SNP[df_SNP$Chr == i,]
        list_snp_chr[[i]]$Chr = gsub("chromosome","chr",list_snp_chr[[i]]$Chr, perl =T)}

return (list_snp_chr);
}

################################################################
# PURPOSE: reading multi files in parallel mode and non parallel
# RETURNS: list of data frames 
################################################################
snpp.readfiles.2Modes <- function(CL,dir_name, reg_exp, N, nCPU, Header, Fill){
  
  files_dir_list = dir(dir_name, full.names = FALSE)    
  indx = grep(reg_exp, files_dir_list, perl = T,  ignore.case = T)
  files_dir_list[indx] = paste(dir_name, files_dir_list[indx], sep ="") 
  
  if(is.null(CL) | nCPU > length(files_dir_list[indx])) {  
      files_list = lapply(files_dir_list[indx], function(my_file){
                      x = try(read.table(my_file, header = Header, fill = Fill))
                      if(class(x) != "try-error") {
                                      return(x)
                       }else{return (-1)}})
    }else{
      files_list = parLapply(CL,files_dir_list[indx], function(my_file){
                      x = try(read.table(my_file, header = Header, fill = Fill))
                      if(class(x) != "try-error") {
                                      return(x)
                       }else{return (-1)}})
      }
  names(files_list) = files_dir_list[indx]
  return (files_list);
}

