rm(list = ls())

args <- commandArgs(trailingOnly=T)
if(length(args) > 0){
  source_dir        = args[1] #"/media/Hasani/SNP_Pipeline/"
  TAS_dir           = args[2] #"/media/Hasani/SNP_Pipeline/"
  output_dir        = args[3] #"/media/Hasani/SNP_Pipeline/"
  job.mode          = args[4]
  nCPU              = args[5]
  joinBlock         = args[6]
  W_basis           = args[7] # either (Pubmed_ID or ACS)
} else{stop (paste("Missing input.", "argument number", length(args)))}

library(plyr)          
library(snow)           # for parallel version of sapply == parSapply
library(ggplot2)

options(warn = 1)
#source_dir        = "/media/Hasani/SNP_Pipeline/scripts/"
#TAS_dir           = "/media/Hasani/SNP_Pipeline/chnge_tas_filter.test/filtered/"
#output_dir        = "/media/Hasani/SNP_Pipeline/chnge_tas_filter.test/filtered/"
#job.mode          = "SOCK"
#nCPU              = 3
#W_basis           = "Pubmed_ID"

print(paste("source_dir", source_dir))
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
clusterEvalQ(CL, source(paste(source_dir,"Utlis.R", sep="")))
clusterEvalQ(CL, source(paste(source_dir,"SNPP_function.R", sep="")))

print ("Reading TAS")
tas_list = snpp.readfiles.2Modes(CL = CL, dir_name = TAS_dir, reg_exp = ".*chr[\\dXYM]*\\.map\\.block$",  N = 9 , nCPU = nCPU, Header = T, Fill=F)

# name both lists by chr
names(tas_list) =unlist(lapply(names(tas_list), function(my_names){
  l = gsub(".+chr","chr", my_names, perl =T,  ignore.case = T)
  l = gsub("\\.map\\.block","", l, perl =T,  ignore.case = T)
  return (l)
}))

print("Found chromosomes:")
print(names(tas_list))
# lapply tas_chr
if(nCPU > length(tas_list)){
  snp_chr_df = lapply(tas_list, function(chr_df){
    f = lapply(chr_df$SNP_Pos, function(pos){
      snp = chr_df[chr_df$SNP_Pos == pos,]
      
      pubID = unique(snp$Pubmed_ID)
      PubID_W = length(pubID)
      
      ACS = unlist(strsplit(as.character(snp$ACS), ","))  
      ACS_W = length(ACS)
      
      df = data.frame(SNP_ID = snp$SNP_ID[1], SNP_Pos = pos, PubID_W =PubID_W, ACS_W =ACS_W, Pubmed_ID= paste(pubID, collapse=","), ACS= paste(ACS, collapse=","), Chr = snp$Chr[1], sPos = snp$sPos[1], ePos = snp$ePos[1], LDP_list = snp$LDP_list[1])
      
      return (df)  })
    
    f = ldply(f, data.frame)
    my_col <- f$SNP_Pos
    my_unique_snp <- f[!(duplicated(my_col)),]
    return(my_unique_snp) })
  
}else{
  print("Paralleling over chromosomes")
  snp_chr_df = parLapply(CL, tas_list, function(chr_df){
    f = lapply(chr_df$SNP_Pos, function(pos){
      snp = chr_df[chr_df$SNP_Pos == pos,]
      
      pubID = unique(snp$Pubmed_ID)
      PubID_W = length(pubID)
      
      ACS = unlist(strsplit(as.character(snp$ACS), ","))  
      ACS_W = length(ACS)
      
      df = data.frame(SNP_ID = snp$SNP_ID[1], SNP_Pos = pos, PubID_W =PubID_W, ACS_W =ACS_W, Pubmed_ID= paste(pubID, collapse=","), ACS= paste(ACS, collapse=","), Chr = snp$Chr[1], sPos = snp$sPos[1], ePos = snp$ePos[1], LDP_list = snp$LDP_list[1])
  
      return (df)  }) 
    
    f = ldply(f, data.frame)
    my_col <- f$SNP_Pos
    my_unique_snp <- f[!(duplicated(my_col)),]
    return(my_unique_snp) })
  }

input_tas_chr= lapply(tas_list, function(l)return(length(unique(l$SNP_Pos))))
names(input_tas_chr) = names(tas_list)

weighted_tas = lapply(snp_chr_df, function(l)return(length(unique(l$SNP_Pos))))
names(weighted_tas) = names(snp_chr_df)

print("Total number of input TAS: ")
print(paste(names(input_tas_chr) , input_tas_chr, sep = " = "))

print("Total number of weighted TAS: ")
print(paste(names(weighted_tas) , weighted_tas, sep = " = "))

my_df1 = t(do.call("rbind", input_tas_chr))
my_df2 = t(do.call("rbind", weighted_tas))

my_df = rbind(my_df1,my_df2)
my_df = transform("file" = c("Input TAS", "Weighted TAS"), my_df)

snpp.saveOutput (my_df, paste(output_dir,"_Read_ME_TAS_weighting_summary.csv", sep=""))

#$TAS_chr1
#SNP_ID   SNP_Pos  Pubmed_ID ACS  Chr  sPos  ePos   LDP_list PubID_W ACS_W

# save the output
if(joinBlock == TRUE){
  print("Warrning: in joinBlock mode, the name of the output file should be given!")
  snp_df = ldply(snp_chr_df, data.frame)
  #Pubmed_ID or ACS
  if(W_basis == "Pubmed_ID"){
    colnames(snp_df)[4] = "SNP_W"
  }else{
    colnames(snp_df)[5] = "SNP_W"
  }
  colnames_out_vect = c("SNP_ID","SNP_Pos","SNP_W", "Chr", "sPos", "ePos", "LDP_list")
  tas_w_df = subset(snp_df, select = colnames_out_vect)
  snpp.saveOutput (tas_w_df, paste(output_dir,".block", sep=""))
  
}else{
  
  output_name_vect = paste(output_dir, names(snp_chr_df), "_W.block", sep="")  
  x = mapply(function(chr_df, names_vect){
    if(nrow(chr_df) > 0){
      if(W_basis == "Pubmed_ID"){
        colnames(chr_df)[3] = "SNP_W"
      }else{
        colnames(chr_df)[4] = "SNP_W"
        }      
        colnames_out_vect = c("SNP_ID","SNP_Pos","SNP_W", "Chr", "sPos", "ePos", "LDP_list")
        tas_w_df = subset(chr_df, select = colnames_out_vect)
        snpp.saveOutput (tas_w_df, names_vect)
    }    
  }, snp_chr_df, output_name_vect)  
}

# Create the list of ACS_w_PubID
tas_df = ldply(tas_list, data.frame)

ACS_vect = unique(unlist(strsplit(as.character(tas_df$ACS), ",")))
ACS_list = vector("list", length(ACS_vect))
names(ACS_list) = ACS_vect

ACS_list = mapply(function (acs, acs_pub){
    hit_df = tas_df[grep(acs, tas_df$ACS, ignore.case=T), ]
    acs_pub = unique(hit_df$Pubmed_ID)
    return(acs_pub)
  }, ACS_vect, ACS_list)

acs_w_pubID = data.frame(ACS = names(ACS_list), W = sapply(ACS_list, length), PubID = sapply(ACS_list, paste, collapse =","), row.names=c(1:length(names(ACS_list))))

snpp.saveOutput (acs_w_pubID, paste(output_dir,"ACS_W_PubID", sep=""))
# =======================================
# plot me
pubIDs = unique(tas_df$Pubmed_ID)
pubIDs_w_list = vector("list", length(pubIDs))
names(pubIDs_w_list) = paste("PubID_", pubIDs, sep="")

pubIDs_w_list = mapply(function (id, pubIDs_list){
  hit_df = tas_df[grep(id, paste("PubID_", tas_df$Pubmed_ID, sep=""), ignore.case=T), ]
  pubIDs_list = nrow(hit_df)
  return(pubIDs_list)
}, names(pubIDs_w_list), pubIDs_w_list)  
 
acs = unique(unlist(strsplit(as.character(tas_df$ACS), ",")))
acs_w_list = vector("list", length(acs))
names(acs_w_list) = acs

acs_w_list = mapply(function (id, acs_list){
  hit_df = tas_df[grep(id, tas_df$ACS, ignore.case=T), ]
  acs_list = nrow(hit_df)
  return(acs_list)
}, names(acs_w_list), acs_w_list)  

FG_Pub_W_plot = ggplot() + geom_point(aes(x = names(pubIDs_w_list), y = pubIDs_w_list), colour = "steelblue") + theme(axis.title.x = element_text(face="bold", colour="#990000"), axis.title.y = element_text(face="bold", colour="#990000"), axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + xlab("Pubmed ID") + ylab("Total number of appearance")

FG_ACS_W_plot = ggplot() + geom_point(aes(x = names(acs_w_list), y = acs_w_list), colour = "steelblue") + theme(axis.title.x = element_text(face="bold", colour="#990000"), axis.title.y = element_text(face="bold", colour="#990000"), axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + xlab("Platform ID") + ylab("Total number of appearance")

save(FG_Pub_W_plot, FG_ACS_W_plot, file= paste(output_dir,"_Foreground_plot.RData", sep=""))
#my_plot = multiplot(g1, g2, cols=1) + ggtitle("Platforms and publications distribution")

#ggsave(paste(output_dir, "ACS_W_PubID_W.pdf"), plot = my_plot, width = 14) 
# save and close everything
stopCluster(CL)
