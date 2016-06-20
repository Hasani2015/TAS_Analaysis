# load libraries and user functions
#source("/data/bioinf/projects/code/projects_misc/Tools/R/UserFunctions/genomeIntervals_userFunctions.R")
#source("/data/bioinf/projects/code/projects_misc/Tools/R/UserFunctions/genomeIntervals_userFunctions_enrichments.R")  # wrapper functions for calculating and plotting the enrichments
rm(list = ls())
# get input
args <- commandArgs(trailingOnly=T)
if(length(args) > 0){
  source_dir        = args[1] #"/media/Hasani/SNP_Pipeline/"
  dir_background    = args[2] #"/media/Hasani/SNP_Pipeline/BG_files_file"
  dir_foreground    = args[3] #"/media/Hasani/SNP_Pipeline/"
  dir_anno          = args[4] #"/media/Hasani/SNP_Pipeline/"
  output_dir        = args[5] #"/media/Hasani/SNP_Pipeline/"
  job.mode          = args[6]
  nCPU              = args[7]
} else{stop (paste("Missing input.", "argument number", length(args)))}

library(plyr)          
library(snow)           # for parallel version of sapply == parSapply
library(ggplot2)
library(genomeIntervals)

options(warn = 1)

print(paste("source_dir", source_dir))
print(paste("dir_background", dir_background)) 
print(paste("dir_foreground", dir_foreground))
print(paste("output_dir", output_dir))
print(paste("dir_anno", dir_anno))
print(paste("job.mode", job.mode))
print(paste("nCPU", nCPU)) 

source(paste(source_dir,"Utlis.R",sep=""))
source(paste(source_dir,"SNPP_function.R",sep=""))
source(paste(source_dir,"genomeIntervals_userFunctions_enrichments_SNPP.R",sep=""))

# initialize parallel version if needed
CL <- NULL
CL <- my.init.parallelRJob.snow(nCPU=nCPU, mode=job.mode, cl=CL)

clusterExport(CL, "source_dir")
clusterEvalQ(CL, library(plyr))   
clusterEvalQ(CL, library(genomeIntervals))
clusterEvalQ(CL, source(paste(source_dir,"Utlis.R", sep="")))
clusterEvalQ(CL, source(paste(source_dir,"SNPP_function.R", sep="")))
clusterEvalQ(CL, source(paste(source_dir,"genomeIntervals_userFunctions_enrichments_SNPP.R", sep="")))

# Read input file
# from: data/testInput.txt
# Read Annotations
print ("Reading ANNO")

  df_Annot_dir = snpp.readFile(dir_anno, F)
  if(nrow(df_Annot_dir)< as.numeric(nCPU)){
	Anno_list = lapply(as.character(df_Annot_dir[,1]), snpp.readFile, H = T)
  }else{Anno_list = parLapply(CL, as.character(df_Annot_dir[,1]), snpp.readFile, H = T)}
  names(Anno_list) = as.character(df_Annot_dir[,1])
  names(Anno_list) = gsub(".*\\/AS\\/","",names(Anno_list))
print(names(Anno_list))

# Read Foreground 
print ("Reading FG")
  df_FG_dir = snpp.readFile(dir_foreground, F)
  files = as.character(df_FG_dir[,1])
  print(files)
  print ("Reading FG list")
  if(nrow(df_FG_dir) < as.numeric(nCPU)){
	Foreground_list = lapply(files, snpp.readFile, H = T)
  }else{ Foreground_list = parLapply(CL, files, snpp.readFile, H = T)}
  
  names(Foreground_list) = as.character(df_FG_dir[,1])
  names(Foreground_list) = gsub(".*\\/Unique_","",names(Foreground_list))
  fg_names_vect = names(Foreground_list)
  print(names(Foreground_list))

  drop ="id"
  Foreground_list = lapply(Foreground_list, function(chr){return (chr[!(names(chr) %in% drop)])})

# Read Background 
  print ("Reading BG")
  df_randBlocks_dir = snpp.readFile(dir_background, F)  
  Background_list = lapply(as.character(df_randBlocks_dir[,1]), snpp.readBlocks, CL = CL, reg_exp = "^Round\\_\\d+.+\\.block$", N = 7, nCPU=nCPU)
  print("DONE")   

  if(exists ("fg_names_vect")){
    # Check Foreground-Background length 
    if(length(Foreground_list) != length(Background_list)){
      stop("overlapGI.R: Foreground and background size mismatch!")
    }else{names(Background_list) = fg_names_vect}	   
  }else{stop("overlapGI.R: error within foreground list (no names found)!")}   

  sampleN = length(Background_list)
  print(names(Background_list))
  print(names(Anno_list))
# Create GenomeInterval objects
#create GI from the Annotation list

  print ("Converting ANNO to GI")
  if(length(Anno_list) < as.numeric(nCPU)){
	Anno_GI_list = lapply(Anno_list, snpp.createGI, ldp = F)
  }else{Anno_GI_list = parLapply(CL, Anno_list, snpp.createGI, ldp = F)}

  print ("DONE!")

  print ("Converting FG to GI")
  #create GI from the Foreground list
  if(length(Foreground_list) < as.numeric(nCPU)){
	Foreground_GI_list = lapply(Foreground_list, snpp.createGI, ldp=T)
  }else{Foreground_GI_list = parLapply(CL, Foreground_list, snpp.createGI, ldp=T)}

  print ("Converting BG to GI")
  #create GI from the Background list
  Background_GI_list = lapply( Background_list, snpp.callDF, CL = CL, ldp=T, nCPU=nCPU)

  print ("DONE!")
#overlap computation
  print ("I'm about to start computing:")
  #my_source = paste(source_dir,"genomeIntervals_userFunctions_enrichments_SNPP.R",sep="")
  oddsRatios <- my.oddsRatiosForAnnotationList(Anno_GI_list, Foreground_GI_list, Background_GI_list, useSnow= FALSE, nCPU=nCPU, nuclOv=TRUE, minOv=0.9, cl=NULL, mode=job.mode, my_source = source_dir) 
  save(oddsRatios, file=paste(as.character(output_dir),"/oddsRatios.RData",sep=""))  
  currentAnno <- names(Anno_list)
  expNames <- names(Background_list)  
  df <- my.testOdds(oddsRatios[["list.ov.orig"]], oddsRatios[["list.ov.sampled"]], oddsRatios[["list.noOv.orig"]],  oddsRatios[["list.noOv.sampled"]],
                          expNames, currentAnno, sampleN)
  
  write.table(df, file=paste(as.character(output_dir),"/oddsRatios.csv",sep=""), quote=FALSE, sep="\t", col.names = TRUE, row.names=FALSE)
  # plot odds ratios and confidence intervals
  my.plotOdds.confidenceInterval(df, pdfFile=paste(as.character(output_dir),"/oddsRatios.confidenceIntervals.pdf",sep=""), doLog2=FALSE)
  my.plotOdds.confidenceInterval(df, pdfFile=paste(as.character(output_dir),"/oddsRatios.confidenceIntervals.log2.pdf",sep=""), doLog2=TRUE)
  # plot relative and absolute observed overlaps
  my.barplot.Odds.mode(df, pdfFile=paste(as.character(output_dir),"/oddsRatios.relativeOverlap.pdf",sep=""), mode="relativeOv", ylab="probes")
  my.barplot.Odds.mode(df, pdfFile=paste(as.character(output_dir),"/oddsRatios.absoluteOverlap.pdf",sep=""), mode="absOv", ylab="probes")
  
# close cluster
#stopCluster(CL)
