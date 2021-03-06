.subtractScores <- function(exp_bigWig_file, input, verbose = TRUE){
    
    chrom_vec <- levels(seqnames(exp_bigWig_file))
    
    if(verbose) message("\t\t Subtracting input")
    
    chrom_input <- as.vector(seqnames(input))
    chrom_exp <- as.vector(seqnames(exp_bigWig_file))
    
    for(chrom in chrom_vec){
        
        scores_input <- score(input[which(chrom_input == 
                                        chrom)])
        scores_exp <- score(exp_bigWig_file[which(chrom_exp 
                                        == chrom)])
        ## If the number of scores differs, correct length
        if(!isTRUE(all.equal(length(scores_input), 
                        length(scores_exp)))){
            
            maxL <- max(length(scores_input), 
                    length(scores_exp))
            
            ## The experiment scores are shorter than the input
            if(length(scores_input) > length(scores_exp))
            {
                scores_input <- scores_input[-which(
                                !(input[which(chrom_input == chrom)] %over% 
                                            exp_bigWig_file[which(chrom_exp 
                                                                    == chrom)]
                                            ))]
            }else{
                length(scores_input) <- maxL
                scores_input[is.na(scores_input)] <- 0
            }
        }
        
        result_scores <- scores_exp - scores_input
        
        score(exp_bigWig_file[which(chrom_exp == chrom)]
        )[which(scores_exp >= 0)] <- result_scores
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
