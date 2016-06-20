
################################################################
# PURPOSE: takes as input an list of genomeIntervals annotations, a list of original and a list of sampled intervals and calculates genome-wide odd rations
# Allows usage of package snow such that a parallel version of sapply speeds up the process
# The lists must be of the following format:
# list.orig     --> "exp1"        --> "1"         --> genomeIntervals object of length nrow(exp1)
#                                 --> "2"
#                                 --> ... 
#                                 --> "sampleN"
#                --> "exp2"
#                --> "exp3"
# annotationList: list of annotations
# list.orig: list containing the original genomeIntervals objects, each entry in list corresponds to one experiment
# list.sampled: list containing the sampled genomeIntervals objects, each entry in list corresponds to one experiment
# useSnow: use snow to call parallel version of sapply == parSapply 
# nCPU: number of CPUs on which parSapply will be distributed
# nuclOv: If overlap should be returned in either number of nucleotides (function: my.ovNucl.genomeIntervals) or number of intervals which overlap to at least minOV% (function: my.ovInterval.genomeIntervals)
#
# Calculate interval_intersection in a subfunction, but not the complete odds ratios. This saves computation time and the number of overlapping nucleotides can be returned, because
# it is needed for the significance test. 
#
# RETURNS: a list of lists as follows:
#	list.odds.orig		list containing for each experiment E a matrix with the Odds for annotation X, i.e. P(A)/(1-P(A))=(ov/N)/(noOv/N)=ov/noOv
#	list.odds.sampled	list containing for each experiment E and sampled intervals S a matrix with the Odds for annotation X, i.e. P(A)/(1-P(A))=(ov/N)/(noOv/N)=ov/noOv
#       list.odds.ratios	list containing for each experiment E and sampled intervals S a matrix with the Odds Ratios, i.e. (ov_orig/noOv_orig)/(ov_sampled/noOv_sampled) 
#				--> for plotting issues only
#       list.ov.orig		list containing for each experiment E a matrix with number of overlapping nucleotides with annotation X	
#       list.ov.sampled		list containing for each experiment E and sampled intervals S a matrix with number of overlapping nucleotides with annotation X
#       list.noOv.orig		list containing for each experiment E a matrix with number of non-overlapping nucleotides with annotation X
#	list.noOv.sampled	list containing for each experiment E and sampled intervals S a matrix with number of non-overlapping nucleotides with annotation X
my.oddsRatiosForAnnotationList <- function(annotationList, list.orig, list.sampled, useSnow=FALSE, nCPU=1, nuclOv=TRUE, minOv=0.9, cl=NULL, mode="SOCK", my_source=""){
################################################################

print (paste("SNOW",useSnow))
print (paste("nCPU",nCPU))
print (paste("mode",mode))
print (paste("my_source",my_source))


	# initiate cluster for parallel computation
	if(useSnow==TRUE){
		library(snow) # load snow library	
		cl <- my.init.parallelRJob.snow(nCPU=nCPU, mode=mode, cl=cl)
		# Load the genomeIntervals package on all cluster nodes
		clusterEvalQ(cl, library(genomeIntervals)) 
		clusterEvalQ(cl, source("/homes/olymp/hoor.al-hasani/Projects/SNP_Pipeline/scripts/genomeIntervals_userFunctions_enrichments_SNPP.R"))
	}

	###
	# init variables
  expNames <- names(list.orig)
  t <- which(expNames %in% names(list.sampled))

  if(length(t) != length(expNames)){ stop("Names in list.orig and list.sampled differ!") }
  list.ov.orig <- list()  	# list containing for each experiment E a matrix with number of overlapping nucleotides with annotation X
	list.ov.sampled <- list()	# list containing for each experiment E and sampled intervals S a matrix with number of overlapping nucleotides with annotation X
	list.noOv.orig <- list()  	# list containing for each experiment E a matrix with number of non-overlapping nucleotides with annotation X
  list.noOv.sampled <- list()	# list containing for each experiment E and sampled intervals S a matrix with number of non-overlapping nucleotides with annotation X
	list.odds.orig <- list()  	# list containing for each experiment E a matrix with the Odds for annotation X, i.e. P(A)/(1-P(A))=(ov/N)/(noOv/N)=ov/noOv
  list.odds.sampled <- list()	# list containing for each experiment E and sampled intervals S a matrix with the Odds for annotation X, i.e. P(A)/(1-P(A))=(ov/N)/(noOv/N)=ov/noOv
  list.odds.ratios <- list()	# list containing for each experiment E and sampled intervals S a matrix with the Odds Ratios --> for plotting issues only

	###
	# loop over all experiments
  for(i in 1:length(expNames)){
		
		  #### overlaps and odds for original intervals
		  if(useSnow==TRUE){
		    list.ov.orig[[i]] <- parSapply(cl, annotationList, my.ovBlock.genomeIntervals, orig = list.orig[[i]])
        print(paste("orig.ov",list.ov.orig[[i]]))
        print("============")
		  }else{
		    list.ov.orig[[i]] <- sapply(annotationList, my.ovBlock.genomeIntervals, orig =list.orig[[i]])	
		    print(paste("orig.ov",list.ov.orig[[i]]))
		    print("============")
      }
	
		  # number of non-overlapping
	  	Norig <- length(list.orig[[i]]@annotation$seq_name)  # total number of blocks in orig      
		  list.noOv.orig[[i]] <- sapply(list.ov.orig[[i]], function(x){ Norig-x } )  
		  # odds 
      list.odds.orig[[i]] <- unlist(list.ov.orig[[i]])/unlist(list.noOv.orig[[i]])
		  ###
	  	# overlaps and odds for sampled intervals
		  list.odds.sampled[[i]] <- list()
      list.odds.ratios[[i]] <- list()
		  list.ov.sampled[[i]] <- list()
		  list.noOv.sampled[[i]] <- list()
		  print("======list.sampled[[i]]======")
      print(length(list.sampled[[i]]))
		  print("============")
      
		  for(j in 1:length(list.sampled[[i]])){  # changed from for(j in 1:sampleN) if the BG size matters for each expr 
		    	print(paste("Orig", i,"vs. Sample", j, sep=" "))
      
		    	if(useSnow==TRUE){
			        list.ov.sampled[[i]][[j]] <- parSapply(cl, annotationList, my.ovBlock.genomeIntervals, orig =list.sampled[[i]][[j]])
			    }else{
			      list.ov.sampled[[i]][[j]] <- sapply(annotationList, my.ovBlock.genomeIntervals, orig =list.sampled[[i]][[j]])	}
			
			# number of non-overlapping
			Nsampled <- length(list.sampled[[i]][[j]]@annotation$seq_name)
      #print(Nsampled)
			list.noOv.sampled[[i]][[j]] <- sapply(list.ov.sampled[[i]][[j]], function(x){ Nsampled-x  })
			# odds
      list.odds.sampled[[i]][[j]] <- unlist(list.ov.sampled[[i]][[j]])/unlist(list.noOv.sampled[[i]][[j]])

			#print("ready")

			# odds ratios
			list.odds.ratios[[i]][[j]] <- list.odds.orig[[i]]/list.odds.sampled[[i]][[j]]
			idx <- which(is.infinite(unlist(list.odds.ratios[[i]][[j]])))
			list.odds.ratios[[i]][[j]][idx] <- NaN
		}
        }

	if(useSnow==TRUE && mode=="SOCK"){
		stopCluster(cl)
	}

  names(list.odds.orig) <- expNames
	names(list.odds.sampled) <- expNames
	names(list.odds.ratios) <- expNames
	names(list.ov.orig) <- expNames
	names(list.ov.sampled) <- expNames
	names(list.noOv.orig) <- expNames
  names(list.noOv.sampled) <- expNames
  return(list(
	list.odds.orig=list.odds.orig, 
	list.odds.sampled=list.odds.sampled, 
	list.odds.ratios=list.odds.ratios, 
	list.ov.orig=list.ov.orig,
	list.ov.sampled=list.ov.sampled, 
	list.noOv.orig=list.noOv.orig, 
  list.noOv.sampled=list.noOv.sampled))
}


################################################################
# PURPOSE: return number of overlapping blocks
##################################################################
my.ovBlock.genomeIntervals <- function(Anno, orig){
  Nov <- 0
  ov <- interval_overlap(orig, Anno)
  if(length(ov)>0){
    #print(paste("my.orig", orig))
    ov_indx_vect = c(1:length(ov))
    Nov_vect = lapply(ov_indx_vect, count.blocks.GI, my_orig= orig, my_Anno = Anno, ov_list = ov)  # lapply over the overlapped blocks list (orig)
    Nov = sum(unlist(Nov_vect))
  }
  return(Nov)
}

##################################################################
# PURPOSE: return 1 if orig_block has LDP within the intersection
# with Anno_block (for loops)
##################################################################
count.blocks.GI1 <- function(vect_ov_indx_crnt_orig_block, my_orig, my_Anno, ov_list){
  block_int = 0
  orig_ldp_vect = 0 
  ov_Anno_block_list = ov_list[[vect_ov_indx_crnt_orig_block]]
  
  if(length(ov_list)>0){
    my_orig = my_orig[vect_ov_indx_crnt_orig_block] 
    for(i in 1:length(ov_Anno_block_list)){  # loop over the annotation blocks which overlapped with the current orig_block
      intrsct_GI = interval_intersection(my_orig, my_Anno[ov_Anno_block_list[i]])
      block_int = data.frame("sPos" = intrsct_GI@.Data[,1], "ePos" = intrsct_GI@.Data[,2])
      orig_ldp_vect = as.numeric(unlist(strsplit(as.character(my_orig$LDP),",")))

      for (i in 1:length(block_int$sPos)){
        orig_ldp = orig_ldp_vect[orig_ldp_vect >= block_int$sPos[i] & orig_ldp_vect <= block_int$sPos[i]]
        if(length(orig_ldp)>0){return (1);}
      }
   }  
  }
  return(0)
}  

##################################################################
# PURPOSE: return 1 if orig_block has LDP within the intersection
# with Anno_block (ldp as GI)
##################################################################
count.blocks.GI <- function(vect_ov_indx_crnt_orig_block, my_orig, my_Anno, ov_list){
  block_int = 0
  orig_ldp_vect = 0 
  ov_Anno_block_list = ov_list[[vect_ov_indx_crnt_orig_block]]
  
  if(length(ov_Anno_block_list)>0){
    my_orig = my_orig[vect_ov_indx_crnt_orig_block] 
    for(i in 1:length(ov_Anno_block_list)){  # loop over the annotation blocks that overlapped with the current orig_block
      intrsct_GI = interval_intersection(my_orig, my_Anno[ov_Anno_block_list[i]])
      
      # Generate GI out of the ldp
      orig_ldp_vect = as.numeric(unlist(strsplit(as.character(my_orig$LDP),",")))
      orig_ldp_df = data.frame(sPos = (orig_ldp_vect)-1, ePos = orig_ldp_vect)
      
      #ldp at block's right border
      LDP_GI = GenomeIntervals(chromosome = rep(my_orig@annotation$seq_name, nrow(orig_ldp_df)), start = orig_ldp_df$sPos, end = orig_ldp_df$ePos, leftOpen = rep(T, nrow(orig_ldp_df)))
      
      my_int_ov = interval_overlap(intrsct_GI, LDP_GI)
      
      my_int_ov_element_length = lapply(my_int_ov, length)
      if(any(my_int_ov_element_length >0)) return(1)
      }
    }  
  return(0)
}

################################################################
# PURPOSE: return number of overlapping probes/segments, which overlap
# to at least minOv%
my.ovInterval.genomeIntervals <- function(X, orig, minOv){
################################################################
        Nov <- 0
        ov <- get.overlap.both(orig, X, minOv)
	ov <- unlist(lapply(ov, length))
	Nov <- length(which(ov>0))
        return(Nov)
}
################# REPLACE THAT - END !!! ########################################

################################################################
# PURPOSE: Test on significance of orig odds compared to expected odds
# In more detail: Fishers exact test is used to compare the number of observed overlapping nucleotides with annotation X with
# number of expected overlap with annotation X. Expected overlap is defined as the mean of observed overlap for all "sampleN" lists of sampled intervals
# INPUT:
# list.ov.orig			list containing for each experiment E a matrix with number of overlapping nucleotides with annotation X
# list.ov.sampled		list containing for each experiment E and sampled intervals S a matrix with number of overlapping nucleotides with annotation X
# list.noOv.orig		list containing for each experiment E a matrix with number of non-overlapping nucleotides with annotation X
# list.noOv.sampled		list containing for each experiment E and sampled intervals S a matrix with number of non-overlapping nucleotides with annotation X
# --> see function my.oddsRatiosForAnnotationList how lists are computed
#
# expNames: names of experiments to identify entries in lists odds.orig and odds.sampled
# currentAnno: names of annotations in an annotationList for which test should be calculated
# sampleN: number of sampled lists of random intervals
# RETURNS:
# dataframe which can be used to plot odds ratio and 95% confidence intervals of fishers exact test by utilizing ggplot2
my.testOdds <- function(list.ov.orig, list.ov.sampled, list.noOv.orig, list.noOv.sampled, expNames, currentAnno, sampleN){
################################################################
	N<-length(expNames)*length(currentAnno)
        v.pvalue <- vector(mode="numeric", length=N)		# vector containing p-values of Fishers test
	v.pvalueStar <- vector(mode="character", length=N)     # vector containing start notation of p-values of Fishers test
        v.ciMin <- vector(mode="numeric", length=N)		# vector containing minimum of 95% confidence interval of odds ratio calculated in Fishers test
        v.ciMax <- vector(mode="numeric", length=N)		# vector containing maximum of 95% confidence interval of odds ratio calculated in Fishers test
        v.oddsRatio <- vector(mode="numeric", length=N)		# vector containing odds ratio of observed odds compared to expected odds: (ov_orig/noOv_orig)/(ov_expected/noOv_expected)
								# vector containing expected odds are estimated as mean of sampleN background lists
	v.anno <- vector(mode="character", length=N)		# vector containing the current annotation for which test is performed
	v.exp <- vector(mode="character", length=N)		# vector containing the current experiment for which test is performed
	v.ov <- vector(mode="numeric", length=N)		# vector containing the observed number of overlapping nucleotides
	v.noOv <- vector(mode="numeric", length=N)              # vector containing the observed number of non-overlapping nucleotides
	v.expectedOv <- vector(mode="numeric", length=N)        # vector containing the expected number of overlapping nucleotides
        v.expectedNoOv <- vector(mode="numeric", length=N)      # vector containing the expected number of non-overlapping nucleotides

	pvalue2starNotation <- function(p){
		if(p < 0.001){ return("***") }
		else if(p < 0.01){ return("**") }
		else if(p < 0.05){ return("*") }
		else{ return(" ") }
	}
print(currentAnno)
print(expNames)
	pos <- 1
	for(e in 1:length(expNames)){
		currentExp <- expNames[e]
    
		for(c in 1:length(currentAnno)){
			aAnno <- currentAnno[c]
      #print(aAnno)
      #print(list.ov.orig[[currentExp]][aAnno])
      
			cont.table <- matrix(nrow=2, ncol=2, dimnames=list(c("overlap", "noOverlap"), c("orig", "expected")))
			cont.table["overlap","orig"] <- as.vector(list.ov.orig[[currentExp]][aAnno])
			cont.table["noOverlap","orig"] <- as.vector(list.noOv.orig[[currentExp]][aAnno])
			sum.ov <- 0
			sum.noOv <- 0
			for(s in 1:sampleN){
				sum.ov <- sum.ov + as.vector(list.ov.sampled[[currentExp]][[s]][aAnno])
				sum.noOv <- sum.noOv + as.vector(list.noOv.sampled[[currentExp]][[s]][aAnno])
			}
			cont.table["overlap","expected"] <- sum.ov/sampleN       
                        cont.table["noOverlap","expected"] <- sum.noOv/sampleN	
			tmp <- fisher.test(cont.table)
		
			v.pvalue[pos] <- tmp$p.value
			v.pvalueStar[pos] <- pvalue2starNotation(v.pvalue[pos])
			v.ciMin[pos] <- tmp$conf.int[1] 
			v.ciMax[pos] <- tmp$conf.int[2]
			v.oddsRatio[pos] <- as.vector(tmp$estimate)
			v.anno[pos] <- aAnno
			v.exp[pos] <- currentExp
			v.ov[pos] <- cont.table["overlap","orig"] 
			v.noOv[pos] <- cont.table["noOverlap","orig"]
			v.expectedOv[pos] <- cont.table["overlap","expected"]
			v.expectedNoOv[pos] <- cont.table["noOverlap","expected"]

			pos <- pos + 1
		}
	}

	df <- data.frame(	pvalue=v.pvalue, ciMin=v.ciMin, ciMax=v.ciMax, oddsRatio=v.oddsRatio,
				ov=v.ov, noOv=v.noOv, expectedOv=v.expectedOv, expectedNoOv=v.expectedNoOv, 
				experiment=v.exp, anno=v.anno, starNotation=v.pvalueStar, stringsAsFactors = FALSE)
	return(df)
}

################################################################
# PURPOSE: Plot odds ratios as points with 95% confidence intervals == result of Fishers exact test 
# df: input dataframe with following columns: 
# data.frame(       pvalue=v.pvalue, ciMin=v.ciMin, ciMax=v.ciMax, oddsRatio=v.oddsRatio,
#                   ov=v.ov, noOv=v.noOv, expectedOv=v.expectedOv, expectedNoOv=v.expectedNoOv,
#                   experiment=v.exp, anno=v.anno)
# see output of my.testOdds 
# Requires ggplot2 library
my.plotOdds.confidenceInterval <- function(df, pdfFile, doLog2=FALSE, width.dodge=0.9){
################################################################
        
	df.tmp <- df

	yLab = "Odds Ratios"
	# plot log2 values
	if(doLog2==TRUE){
		df.tmp$ciMin <- log2(df.tmp$ciMin)
		df.tmp$ciMax <- log2(df.tmp$ciMax)
		df.tmp$oddsRatio <- log2(df.tmp$oddsRatio)
		yLab="Log2(Odds Ratios)"
	}

	# in case of inifite confidence intervals (either due to log2 transformation or no observed overlap) --> do not plot it
	idx <- unique(which(is.infinite(df.tmp$ciMin) | is.infinite(df.tmp$ciMax) | is.infinite(df.tmp$oddsRatio)))
	if(length(idx)>0){
		df.tmp[idx,]$ciMin <- NaN
		df.tmp[idx,]$ciMax <- NaN
		df.tmp[idx,]$oddsRatio <- NaN
		df.tmp[idx,]$starNotation <- "NaN"
	}
	
        #limits <- aes(ymin = df.tmp$ciMin, ymax=df.tmp$ciMax)  # this does not work, because: gplot will look for the variables in aes in either the global environment (when the dataframe is specifically added as df$... ) 
								# or within the mentioned environment.
        limits <- aes(ymin=ciMin, ymax=ciMax)
	p <- ggplot(df.tmp, aes(colour=experiment, y=oddsRatio, x=anno)) #, group=experiment))
        p <- p + geom_point(position=position_dodge(width=width.dodge), shape=15, size=2) 
	p <- p + geom_errorbar(position=position_dodge(width=width.dodge), limits, width=0.6)
        p <- p + scale_color_manual(values="steelblue")
        p <- p + labs(x="", y=yLab)
	p <- p + geom_text(aes(y=oddsRatio, label=starNotation), hjust=-0.3, position=position_dodge(width=width.dodge))
        p <- p + opts( axis.text.x=theme_text(angle=90, size = 12),
                        axis.text.y=theme_text(size = 12),
                        #legend.text=theme_text(size = 12), 
                        axis.title.y=theme_text(angle=90, size = 15)
                        #legend.title=theme_text(size = 0),
                        #legend.position="none"
		)
        #p <- p + facet_wrap(~experiment, ncol=1, nrow=3, scales="free_y")

	# set vertical lines such that groups are visually separated
	nAnno <- length(unique(df.tmp$anno))
	offsetVLine <- 0.25 
	p <- p + geom_vline(xintercept = 1+offsetVLine:nAnno+offsetVLine, colour="grey", size=0.5)
#	p <- p + geom_vline(xintercept = 1-offsetVLine:nAnno-offsetVLine, colour="grey", size=0.5)
	# set horizontal line denoting "no enrichment" and "no de-enrichment"
	p <- p + geom_abline(intercept=0, slope=0, colour="black", size=0.5)
 
	#print(p)
	ggsave(pdfFile, plot = p)
        svgFile<-gsub(".pdf", ".svg", pdfFile)
	#print(p)
	ggsave(svgFile, plot = p)        
}

################################################################
# PURPOSE: Barplot with different values according to chosen mode 
# df: input dataframe with following columns: 
# data.frame(       pvalue=v.pvalue, ciMin=v.ciMin, ciMax=v.ciMax, oddsRatio=v.oddsRatio,
#                   ov=v.ov, noOv=v.noOv, expectedOv=v.expectedOv, expectedNoOv=v.expectedNoOv,
#                   experiment=v.exp, anno=v.anno)
# modes: 
# relativeOv    	barplot of observed relative overlaps
# absOv         	barplot of observed absolute overlaps
# relativeExpOv		barplot of expected relative overlaps
# absExpOv		barplot of expected absolute overlaps
#
# see output of my.testOdds 
# Requires ggplot2 library
my.barplot.Odds.mode <- function(df, pdfFile="", mode="relativeOv", ylab="nucleotides"){
################################################################
        
	if(!mode %in% c("relativeOv", "absOv", "relativeExpOv", "absExpOv")){
		stop(paste("Unknown mode, should be one of: ", c("relativeOv", "absOv", "relativeExpOv", "absExpOv"), sep=" "))
	}

	df[, "val"] <- vector(mode="numeric", length=nrow(df))
	ylabString=""
	if(mode=="relativeOv"){
		df[,"val"] <- df[,"ov"]/(df[,"ov"]+df[,"noOv"])
		ylabString=paste("Fraction of overlapping ", ylab, sep="")
	}else if(mode=="absOv"){
		df[,"val"] <- df[,"ov"]
		ylabString=paste("Number of overlapping ", ylab, sep="")
	}else if(mode=="relativeExpOv"){
                df[,"val"] <- df[,"expectedOv"]/(df[,"expectedOv"]+df[,"expectedNoOv"])
		ylabString=paste("Expected fraction of overlapping ", ylab, sep="")
        }else if(mode=="absExpOv"){
                df[,"val"] <- df[,"expectedOv"]
		ylabString=paste("Expected number of overlapping ", ylab, sep="")
        }else{
		stop(paste("Unknown mode, should be one of: ", c("relativeOv", "absOv", "relativeExpOv", "absExpOv"), sep=" "))
	}

        p = ggplot()
        p = p + geom_bar(aes(anno,val, fill = experiment), data = df)
        #p <- ggplot(df, aes(anno, val, fill="steelblue"))
        #p <- p + geom_bar(position="dodge", stat="identity")
        #p <- p + scale_fill_manual(values=factor("steelblue"))
        p <- p + labs(x="", y=ylabString)
        p <- p + opts( axis.text.x=theme_text(angle=90, size = 12),
                       axis.text.y=theme_text(size = 12),
                       axis.title.y=theme_text(angle=90, size = 15))
               
	# PDF file
  ggsave(pdfFile, plot = p)
        #pdf(file=pdfFile)
        #print(p)
       # dev.off()
        # SVG file
        svgFile<-gsub(".pdf", ".svg", pdfFile)
        #svg(file=svgFile)
       # print(p)
      #  dev.off()
ggsave(svgFile, plot = p)

}


################################################################
# PURPOSE: keep only those entries with valid seq_names, i.e. from chr1 to chrM, otherwise
# assignChrRanges fails after calling interval_complement
my.validSeqNames <- function(x){
################################################################
	validSeqNames <- as.vector(chrRanges$seq_name)
	x <- x[which(seq_name(x) %in% validSeqNames),]
	x@annotation$seq_name <- as.vector(x@annotation$seq_name)
	x@annotation$seq_name <- as.factor(x@annotation$seq_name)
	return(x)
}


