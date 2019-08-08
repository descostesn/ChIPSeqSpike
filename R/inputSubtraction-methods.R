.subtractScores <- function(exp_bigWig_file, input, verbose = TRUE){
    
    chrom_vec <- levels(seqnames(exp_bigWig_file))
    
    if(verbose) message("\t\t Subtracting input")
    
    chrom_input <- as.vector(seqnames(input))
    chrom_exp <- as.vector(seqnames(exp_bigWig_file))
    
    for(chrom in chrom_vec){
        
		if(verbose) message("\t\t\t processing chrom: ", chrom)
		
		gr_exp_chrom <- exp_bigWig_file[which(chrom_exp == chrom)]
		gr_input_chrom <- input[which(chrom_input == chrom)]
		
		all_scores_input <- score()
		all_scores_exp <- score(gr_exp_chrom)
		all_start_input <- start(gr_input_chrom)
		all_start_exp <- start(gr_exp_chrom)
		all_end_input <- end(gr_input_chrom)
		all_end_exp <- end(gr_exp_chrom)
		all_strand_input <- as.character(strand(gr_input_chrom))
		all_strand_exp <- as.character(strand(gr_exp_chrom))
		all_score_input <- score(gr_input_chrom)
		all_score_exp <- score(gr_exp_chrom)
		
		
		## From the overlap below, several ranges of input can be overlapping 
		## an experiment range. The logic is to subtract each input score to
		## the same exp score and build a new grange for exp including the
		## sub-intervals
		
		result_overlap <- findOverlaps(gr_exp_chrom, gr_input_chrom)
		
		input_index_list <- split(subjectHits(result_overlap), 
				factor(queryHits(result_overlap)))
		
		exp_index_list <- as.list(unique(queryHits(result_overlap)))
		
		if(!isTRUE(all.equal(length(input_index_list), length(exp_index_list))))
			stop("This should not happen during input subtraction, contact ",
					"the developper.")
		
		length_vec <- unlist(lapply(input_index_list, length))
		more_than_one_nb <- sum(length_vec[which(length_vec > 1)])
		new_exp_granges <- vector(mode="list", length = more_than_one_nb)
		
		!! To optimize !!!!!!!!!!!!!!
		
		start_time <- Sys.time()
		#for(i in seq_len(length(exp_index_list))){
		for(i in seq_len(10)){	
			input_index_vec <- input_index_list[[i]]
			exp_index <- exp_index_list[[i]]
			input_index_length <- length(input_index_vec)
			
			if(input_index_length > 1){
				
				new_exp_scores <- all_scores_exp[exp_index]-
						all_scores_exp[input_index_vec]
				
				new_exp_start <- pmax(rep(all_start_exp[exp_index], 
								input_index_length), 
						all_start_input[input_index_vec])
				
				new_exp_end <- pmin(rep(all_end_exp[exp_index],
								input_index_length),
						all_end_input[input_index_vec])
				
				new_exp_chrom <- rep(chrom, input_index_length)
				new_exp_strand <- rep(all_strand_exp[exp_index], 
						input_index_length)
				
				new_exp_granges[[exp_index]] <- GRanges(
						seqnames = new_exp_chrom,
						ranges = IRanges(start = new_exp_start,
								end = new_exp_end),
						strand = new_exp_strand,
						score = new_exp_scores)
			}else{
				new_exp_granges[[exp_index]] <- gr_exp_chrom[exp_index]
				score(gr_exp_chrom[exp_index]) <- all_score_exp[exp_index] -
						all_score_input[input_index_vec]
			}
		}
		end_time <- Sys.time()
		end_time - start_time
		
#		scores_input <- unlist(lapply(input_index_list, 
#						function(index_vec) 
#							return(max(all_scores_input[index_vec]))))
#		scores_exp <- score(exp_bigWig_file[which(chrom_exp == chrom)])
#		
#		result_scores <- scores_exp - scores_input
#		score(exp_bigWig_file[which(chrom_exp == chrom)]) <- result_scores
		
    }
    
    return(exp_bigWig_file)
}


setMethod(
        
        f = "inputSubtraction",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, verbose = TRUE){
            
            if(.Platform$OS.type != 'windows') {
                if(!length(grep("RPM", getBigWigFile(theObject))))
                    stop("RPM normalization must be performed before ",
                            "subtracting the input")
                
                input_bigWig_file <- getLoadedData(theObject)
                
                if(verbose) message("Subtracting input to experiment")
                
                experimentList(theObject) <- lapply(
                        getExperimentList(theObject), function(exp, input){
                            
                            if(verbose)
                                message("\t Processing ", getExpName(exp))
                            
                            exp_bigWig_file <- getLoadedData(exp)
                            
                            seqnanamesBW <- seqnames(exp_bigWig_file)
                            if(!isTRUE(all.equal(levels(seqnanamesBW),
                                            levels(seqnames(input)))))
                                stop("Chromosomes differ between input and ",
                                        "experiment")
                            
                            if(!isTRUE(all.equal(
                                                    as.character(
                                                            seqnames(input)[1]
                                    ),
                                                    as.character(seqnames(
                                                        exp_bigWig_file)[1])))
                                    || !isTRUE(all.equal(as.numeric(
                                                            start(input)[1]),
                                                    as.numeric(
                                                            start(
                                                        exp_bigWig_file)[1])))
                                    || !isTRUE(all.equal(
                                                    as.numeric(end(input)[1]), 
                                                    as.numeric(end(
                                                        exp_bigWig_file)[1]))))
                    stop("The binning coordinates between the input DNA",
                            " and the experiment are different. The Bigwig ",
                            "files were probably generated using a different",
                            " procedure.")
                            
                            if(!isTRUE(all.equal(as.numeric(end(input)[2])-
                                                as.numeric(start(input)[2]),
                                        as.numeric(end(exp_bigWig_file)[2])
                        -as.numeric(start(exp_bigWig_file)[2]))))
                                stop("The bin size is different between the",
                                        " input DNA and the experiment. The ",
                                        "Bigwig files were probably generated",
                                        " using a different procedure.")
                            
                            result <- .subtractScores(exp_bigWig_file, 
                                    input, verbose)
                            
                            loadedData(exp) <- result
                            bigWigFile(exp) <- .modifyBigWigName(exp, "BGSub")
                            
                            return(exp)
                            
                        }, input_bigWig_file)
                
                return(theObject)
            }else{
                stop("As of rtracklayer >= 1.37.6, BigWig is not supported",
                        " on Windows.")
            }
        })


setMethod(
        
        f = "inputSubtraction",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, verbose = TRUE){
            
            if(.Platform$OS.type != 'windows') {
                if(!length(grep("RPM", getBigWigFile(theObject))))
                    stop("RPM normalization must be performed before ",
                            "subtracting the input")
                
                if(verbose) message("Reading the input file")
                
                input_file_path <- getBigWigFile(theObject)
                
                input_bigWig_file <- import(input_file_path, format="BigWig")
                
                if(verbose) message("Subtracting input to experiment")
                
                experimentList(theObject) <- lapply(
                        getExperimentList(theObject), function(exp, input){
                            
                            if(verbose){
                                
                                message("\t Processing ", getExpName(exp))
                                message("\t\t Reading bigWig file")
                            }
                            
                            exp_bigWig_file <- import(getBigWigFile(exp), 
                                    format="BigWig")
                            
                            if(!isTRUE(all.equal(levels(
                                                    seqnames(exp_bigWig_file)),
                                            levels(seqnames(input)))))
                                stop("Chromosomes differ between input and ",
                                        "experiment")
                            
                            if(!isTRUE(all.equal(as.character(
                                                            seqnames(input)[1]
                                    ),
                                    as.character(seqnames(
                                                    exp_bigWig_file)[1])))
                                    || !isTRUE(all.equal(as.numeric(
                                                            start(input)[1]),
                                                    as.numeric(
                                                            start(
                                                        exp_bigWig_file)[1])))
                                    || !isTRUE(all.equal(as.numeric(
                                                            end(input)[1]), 
                                                    as.numeric(end(
                                                        exp_bigWig_file)[1]))))
                    stop("The binning coordinates between the input DNA",
                            " and the experiment are different. The ",
                            "Bigwig files were probably generated using ",
                            "a different procedure.")
                            
                            if(!isTRUE(all.equal(
                                            as.numeric(
                                                    end(input)[2])-
                                            as.numeric(start(input)[2]),
                                    as.numeric(end(exp_bigWig_file)[2])-
                                            as.numeric(start(
                                                        exp_bigWig_file)[2]))))
                                stop("The bin size is different between the",
                                        " input DNA and the experiment. The ",
                                        "Bigwig files were probably generated",
                                        " using a different procedure.")
                            
                            result <- .subtractScores(exp_bigWig_file, input, 
                                    verbose)
                            
                            
                            output_bigWig <- .modifyBigWigName(exp, "BGSub",
                                    NULL)
                            export(result, con = output_bigWig, 
                                    format="BigWig")
                            
                            bigWigFile(exp) <- output_bigWig
                            return(exp)
                            
                        }, input_bigWig_file)
                
                return(theObject)
            }else{
                stop("As of rtracklayer >= 1.37.6, BigWig is not supported on",
                        "Windows.")
            }
        }
)


.loadInputSubtractionList <- function(theObject, verbose){
    
    datasetList(theObject) <- lapply(getDatasetList(theObject), 
            function(object){
                
                return(inputSubtraction(object, verbose))
            })
    
    return(theObject)
}


setMethod(
        
        f = "inputSubtraction",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject, verbose = TRUE){
            
            if(.Platform$OS.type != 'windows') {
                .loadInputSubtractionList(theObject, verbose)
            }else{
                stop("As of rtracklayer >= 1.37.6, BigWig is not supported on",
                        " Windows.")
            }
        }
)


setMethod(
        
        f = "inputSubtraction",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject, verbose = TRUE){
            
            if(.Platform$OS.type != 'windows') {
                .loadInputSubtractionList(theObject, verbose)
            }else{
                stop("As of rtracklayer >= 1.37.6, BigWig is not supported",
                        " on Windows.")
            }
        }
)
