#' @title Identify crossover, non-crossover, and gene conversion tracts
#' 
#' @description Infers different types of tracts (crossover, non-crossover, and gene conversion) along a 
#' chromosome given genotyped yeast tetrads or simulated data
#' 
#' @param data a \code{states.matrix} or \code{forward.backward} object inherited from either 
#' \code{recombine_to_tetrad_states} or \code{estimate_anc_fwd_back} 
#' 
#' @param tetrad (Optional) A numeric or character specifying the tetrad ID
#' 
#' @param chr (Optional) A numeric or character specifying the chromosome ID if class \code{tetrad.states}. 
#' If \code{data} is of class \code{forward.backward} than \code{chr} will be provided.
#' 
#' @details (to do)
#' 
#' @return A data.frame containing the following columns:
#' \describe{
#'      \item{tetrad}{The tetrad ID (default == 1)}
#'      \item{chr}{The chromosome ID (default == "I")}
#'      \item{type}{The type of inferred tract:} 
#'      \describe{
#'          \item{1. 2_2}{a tract with a 2:2 ratio}
#'          \item{2. GC_tel}{a gene conversion tract located at a chromosomal end}
#'          \item{3. GC_internal}{a gene conversion tract with at least one flanking tract also a gene conversion tract}
#'          \item{4. COnoGC}{indicates where a crossover event occurred but without a detectable gene conversion tract}
#'          \item{5. COyesGC}{a gene conversion tract associated with a crossover event}
#'          \item{6. NCO}{a gene conversion tract without an associated crossover event}
#'      }
#'      \item{start_snp}{numeric value of the first snp in a given tract}
#'      \item{end_snp}{numeric value of the last snp in a given tract}
#'      \item{extent}{numeric value of the range of snps (inclusive) in a given tract}
#'      \item{bias}{integer specifying the gene conversion bias where:}
#'      \describe{
#'          \item{a. 0}{a 4:0 segregation (all of parent 0 type)}
#'          \item{b. 1}{a 3:1 segregation}
#'          \item{c. 2}{a 2:2 no bias found}
#'          \item{d. 3}{a 1:3 segregation}
#'          \item{e. 4}{a 0:4 segregation (all of parent 1 type)}
#'          \item{f. NA}{(COnoGC only) bias information is not applicable}  
#'      } 
#' }
#' 
#' @seealso \code{\link{recombine_index}}, \code{\link{recombine}}, \code{\link{recombine_to_tetrad_states}}
#' 
#' @author Tyler D. Hether
#' 
#' @export infer_tracts
#' 
#' @examples
# Example 1: 1 tetrad
set.seed(1) # For reproducability
l <- 1000 # number of loci to simulate
rec <- 0.01 # recombination rate between each snp

r <- recombine_index(rec, 1:l) # recombination rate between each snp (vector form)
p_a <- .999 # probability of correct sequencing assignment
p <- make_parents(floor(seq(from=1, to=1e5, length.out=l))) # make the parent
recomb_sim <- recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert=1, length.conversion=10) # recombine parents
states <- recombine_to_tetrad_states(tetrad_data=recomb_sim) # convert to tetrad.states object

threshold_size <- 1e4

df <- ddply(states, .(Tetrad, Chr), infer_tracts, threshold_size=1e2)
hist(dplyr::filter(df, type=="COyesGC" | type=="COnoGC" | type=="NCO")$start_snp,
 breaks=200, xlab="snp", main="start position of recombination point")
#
#' # Example 2: 100 simulated tetrads with a recombination hotspot
#' set.seed(1) # For reproducability
#' rec <- c(rep(0.001, 99), 0.4, rep(0.001, 99))
#' res <- sim_tetrad(n.tetrads=250, l=200, rec=rec, p.assign=0.999, 
#'    mu.rate=0, f.cross=0.8, f.convert=0, length.conversion=0, coverage=0.5)
#' snp.dat <- tetrad_to_df(res)
#' # Get the inferred states
#' states1 <- ddply(snp.dat, .(Tetrad, Spore, Chr), function(x){
#'     est_fwd_back(snp.dat=x, p.assign=0.999, p.trans=0.01)
#'     })
#' states2 <- fwd_back_to_tetrad_states(fb_data=states1)
#' df1 <- ddply(states2, .(Tetrad, Chr), infer_tracts)
#' 
#' hist(dplyr::filter(df1, type=="COyesGC" | type=="COnoGC" | type=="NCO")$start_snp,
#'  breaks=200, xlab="snp", main="start position of recombination point")

infer_tracts <- function(data, threshold_size=1e4){
    #! require(reshape2)
   
    if(!("Tetrad" %in% colnames(data))){
        stop("No Tetrad column detected")
         # inferred_states_data$Tetrad <- as.numeric(unlist(lapply(strsplit(inferred_states_data$Ind, split="_"), "[", 1L)))
    } 
    # if(!("Spore" %in% colnames(inferred_states_data))){
    #     warning("No Spore column detected, trying to infer it from Ind column")
    #      inferred_states_data$Spore <- as.numeric(unlist(lapply(strsplit(inferred_states_data$Ind, split="_"), "[", 2L)))
    # } 

    #! data <- dcast(inferred_states_data, Tetrad + Chr + Snp ~ Spore, value.var="states_inferred", fun.aggregate = mean)

    # data is a data.frame with the following columns
    colnames(data) <- c("Tetrad", "Chr", "Snp", "one", "two", "three", "four")
    # Order the data
    data <- data[with(data, order(Tetrad, Chr, Snp)),]
    # filter out any NaNs
    data <- data[complete.cases(data), ]


    # check_dat(data)
    # # Input is a unique combination of tetrad id and chr id.
    # if(length(unique(data[,1]))!=1){
    #     stop("only a unique tetrad id should be passed to infer_tracts")
    # }
    # if(length(unique(data[,2]))!=1){
    #     stop("only a unique chromosome id should be passed to infer_tracts")
    # }    

    OUT <- ddply(data, .(Tetrad, Chr), function(tetrad_states, ...){

        # tetrad_states <- data
        tetrad <- tetrad_states$Tetrad[1]
        chr <- tetrad_states$Chr[1]

        # Get the unique candidate 'regions':
        regions <- unique_regions(tetrad_states)

        # Calculate the biases, if present:
        biases <- check_GC_bias(tetrad_states)

        # Get the actual segregation pattern:
        text <- get_text(tetrad_states)

        # Combine datasets: 
        CO_summary <- data.frame(snp=regions$tetrad_states[,'Snp'], type=regions$type, GCbias=biases, text=text)


        # unique(CO_summary$type)
        # Go through the chromosome and call each region as a specific 'tract':
        out <-   do.call(rbind, lapply(unique(CO_summary$type), function(res, ...){
            # print(res)
            type="NULL"
            BIAS="NULL"
            # print(res)
            res_range <- range(CO_summary$snp[which(CO_summary$type==res)])
            extent <- 1+(res_range[2] - res_range[1])
            # Get the index position for this region.
            pos <- which(unique(CO_summary$type)==res)
            
            # Is this region smaller than the threshold size?
            if(extent<threshold_size){
                smaller_than_threshold <- TRUE
            } else {
                smaller_than_threshold <- FALSE
            }

            # If the region is greater than the threshold *and* has a 2:2 pattern, 
            # call it a "2_2" region. Note that COnoGC regions will be checked
            # for at the end. Note that 2:2 tracts that are smaller than the 
            # threshold size will be considered non-2:2 (GC) tracts for the rest of the
            # algorithm.
            if(!smaller_than_threshold & res>0){
                type <- "2_2"
                # Return the degree of bias
                BIAS <- 2
            }

            # if the region is smaller than the threshold *or* it was previously considered
            # a non 2:2 tract, do the following: 
            if(res < 0 | smaller_than_threshold){
                # Is this 'non 2:2' (GC) region on a chromosomal end? 
                # If so, classify it as 'GC_tel'. If not, temporarily call it 'internal'.
                if(res==unique(CO_summary$type)[1] | res==rev(unique(CO_summary$type))[1]){
                    # This block is executed if the region *is* on the ends:
                    # internal <- FALSE
                    type <- "GC_tel"
                } else {
                    # internal <- TRUE
                    # If this 'non 2:2' (GC) region is not internal, identify whether it
                    # is part of a NCO event or part of a CO event. To do this, first identify
                    # the left and right flanking regions, their extent, and whether they have 
                    # been classified as "GC" or "nonGC". "GC" is anything region that is 
                    # non 2:2 or is 2:2 but smaller than the threshold_size.
                    left <- get_flanking(res_range=res_range, direction="left", df=CO_summary, threshold_size=threshold_size)
                    right <- get_flanking(res_range=res_range, direction="right", df=CO_summary, threshold_size=threshold_size)                
                    
                    # Second, if either of the flanking regions are GC regions (non 2:2 or 2:2 <
                    # threshold_size) then temporarily call it a GC_internal (to be reclassified
                    # later). If both of the flanking regions are "nonGC" then determine NCO (identical
                    # flanking pattern) or COyesGC (i.e., different flanking pattern)
                    if(left$type2=="GC" | right$type2=="GC"){
                        type <- "GC_internal"
                    } else {
                        if(as.character(left$text) == as.character(right$text)){
                            type <- "NCO"
                        } else {
                            type <- "COyesGC"
                        }
                    }
                }
                # Return the degree of BIAS
                BIAS <- CO_summary[CO_summary$type==unique(CO_summary$type)[pos],'GCbias'][1]
            }

            # Store results as a data.frame
            return(data.frame(res, tetrad, chr, type, res_range[1],res_range[2], extent,BIAS,
                pattern=as.character(CO_summary[CO_summary$type==unique(CO_summary$type)[pos],'text'][1])))
        }))
        
        # Return output of inferred tracts along the chromosome 
        # out <- out[,-1] # Housekeeping
        colnames(out) <- c("region", "tetrad", "chr", "type", "start_snp", "end_snp", "extent", "bias", "pattern")
        class(out) <- c("data.frame", "inferred.tracts")

        # Reclassifying temporary tracts:

        # This identifies the approx location of COs without a detected gene conversion tract
        # print("Screening for COs without GCs...")
        out2 <- infer_COnoGC_tracts(out)

        # Now merge GC_internal tracts, if appropriate:
        # print("Merging internal GC tracts...")
        
        i=0
        out3 <- NULL
        location_of_2_2s <- which(out2$type[1:dim(out2)[1]]=="2_2")
        while(i<dim(out2)[1]){
            i <- i + 1
            # print(i)
            # If region i is GC of some kind...
            if(out2$type[i]!="2_2"){
                # ...and there is another 2:2 tract 3' along the chromosome...
                if(length(which(location_of_2_2s>i))>0){
                    # ...find the last_tract_to_merge index...
                    last_tract_to_merge <- location_of_2_2s[location_of_2_2s>i][1]-1
                    # ...and merge the ith row with the last_tract_to_merge-th row. If 
                    # one of the GC tracts is labeled GC_tel, convert the whole tract
                    # to a GC_tel. 
                    if("GC_tel" %in% out2$type[i:last_tract_to_merge]){
                         out3 <- rbind(out3, data.frame(region=out2$region[i], 
                                tetrad=out2$tetrad[i], 
                                chr=out2$chr[i],
                                type="GC_tel",
                                start_snp=out2$start_snp[i], 
                                end_snp=out2$end_snp[last_tract_to_merge])) 
                         # And reset the counter
                         i <- last_tract_to_merge          
                    } else {
                        # If, however, all GC tracts in this group were internal then
                        # determine if it was a NCO or a COyesGC
                        if(as.character(out2$pattern[i-1])==as.character(out2$pattern[last_tract_to_merge+1])){
                            # if the flanking 2:2 patterns (>the threshold) are the same
                            # then call the whole, merged region an NCO...
                             out3 <- rbind(out3, data.frame(region=out2$region[i], 
                                    tetrad=out2$tetrad[i], 
                                    chr=out2$chr[i],
                                    type="NCO",
                                    start_snp=out2$start_snp[i], 
                                    end_snp=out2$end_snp[last_tract_to_merge])) 
                        } else {
                            # ...otherwise, call the whole region a COyesGC
                             out3 <- rbind(out3, data.frame(region=out2$region[i], 
                                    tetrad=out2$tetrad[i], 
                                    chr=out2$chr[i],
                                    type="COyesGC",
                                    start_snp=out2$start_snp[i], 
                                    end_snp=out2$end_snp[last_tract_to_merge])) 
                        }
                        # And reset the counter
                        i <- last_tract_to_merge   
                    }
                }  else {
                    # ...but if there wasn't another 2:2 tract along the chromosome
                    # then the remainder of the regions are of type GC_tel
                    out3 <- rbind(out3, data.frame(region=out2$region[i], 
                           tetrad=out2$tetrad[i], 
                           chr=out2$chr[i],
                           type="GC_tel",
                           start_snp=out2$start_snp[i], 
                           end_snp=out2$end_snp[dim(out2)[1]])) 
                    # And reset the counter
                    i <- dim(out2)[1]  
                }
            } else {
               # If the region is a 2:2 pattern then just return the row
                out3 <- rbind(out3, data.frame(region=out2$region[i], 
                       tetrad=out2$tetrad[i], 
                       chr=out2$chr[i],
                       type="2_2",
                       start_snp=out2$start_snp[i], 
                       end_snp=out2$end_snp[i]))                
            }
        }

        out3$extent <- (out3$end_snp - out3$start_snp) + 1

        return(out3)
    })
   
    # Housekeeping: 
    OUT <- OUT[,-c(3,4,5)]
    return(OUT)
}

#' @title Identify recombination points from state sequences
#' 
#' @description This is a simple function that takes a vector of parental states
#' along a chromosome and identifies where states have changed. Note: In this method, 
#' a single base pair mutation is indistinguishable from a double recombination event. 
#' 
#' @param state.vector a vector of state values (0 or 1)
#'
#' @return a vector of length \code{state.vector} stating whether a recombination 
#' event occured (1) or not (0).
#' 
#' @seealso \code{\link{sim_en_masse}}
#'
#' @importFrom dplyr filter
#' 
#' @author Tyler D. Hether 
#' 
#' @export id_recombination_events
#' 
#' @examples
#' # Example 1: simple
#' # A recombination occurred between snp 3 and 4 and between 8 and 9.
#' statepath <- c(0,0,0,1,1,1,1,1,0,0,0)
#' id_recombination_events(state.vector=statepath)
#' which(id_recombination_events(state.vector=statepath)==1)
#' #
#' # Example 2: complex
#' set.seed(1) # For reproducability
#' # simulate a recombination hotspot between the 99th and 100th snp
#' rec <- c(rep(0.001, 99), 0.4, rep(0.001, 99))
#' # simulate 500 spores en masse
#' n.spores <- 500 
#' spores <- sim_en_masse(n.spores=n.spores, l=200, rec=rec, 
#'  p.assign=.999, mu.rate=0.001, f.cross=0.5, 
#'     f.convert=0.5, length.conversion=10, coverage=1)
#' # Convert to dataframe
#' snp.dat <- en_masse_to_df(spores)
#' # Infer states
#' states1 <- ddply(snp.dat, .(Tetrad, Spore, Chr), function(x){
#'     est_fwd_back(snp.dat=x, p.assign=0.999, p.trans=mean(rec))
#'     })
#' # ddply through each spore to find recombination points (rpts)
#' df <- ddply(states1, .(Spore), function(x){
#'     rpts <- which(id_recombination_events(x$states_inferred)==1)
#'     npts <- length(rpts)
#'     return(data.frame(rpts=rpts))
#'     })
#' # Plot
#' hist(df$rpts, breaks=200, xlab="snp", 
#'     main="recombination frequency")

id_recombination_events <- function(snps_genotypes_df){ 
        # snps is a vector containing the snp locations 
        # genotypes is a vector of corresponding genotypes
        # colnames(snps_genotypes_df) <- c("snp", "genotype")
        # snps_genotypes_df <- snps_genotypes_df[with(snps_genotypes_df, order(snp)),]


        dims <- dim(snps_genotypes_df)[1]
        if(dims<=1){
            stop("At least two snps needed to infer recombination points.")
        }
        no0yes1 <- sapply(2:dims, function(i){
            if(snps_genotypes_df[i,2]!=snps_genotypes_df[i-1,2]){
                return(1)
            } else 
            {
                return(0)
            }
            })
        midpoints <- sapply(2:dims, function(i){
            return(mean(c(snps_genotypes_df[i,1], snps_genotypes_df[i-1,1])))
        })
    
        return(data.frame(midpoints=midpoints, no0yes1=no0yes1))

}

# Minor functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Checks to make sure input data frame will work for infer_tracts function:
check_dat <- function(dataframe){
    if(dim(dataframe)[2]!=7){
        stop("input data needs 7 columns c('Tetrad', 'Chr', 'Snp', 'p0', 'p1')")
    }
    if(!is.numeric(dataframe[,1]) & !is.integer(dataframe[,1]) & !is.character(dataframe[,1])){
        stop("Tetrad ids need to be numeric, integers or characters")
    }
    if(!is.numeric(dataframe[,2]) & !is.integer(dataframe[,2]) & !is.character(dataframe[,2])){
        stop("Chromosome ids need to be numeric, integers or characters")
    }
    if(!is.numeric(dataframe[,3]) & !is.integer(dataframe[,3])){
        stop("snp ids need to be numeric or integers")
    }
    # Check that snp ids are unique and sorted for a given tetrad:chr combo:
    sort_unique <- ddply(dataframe, .(Tetrad, Chr), function(x){
            snps <- x$Snp
            SORT <- FALSE
            UNIQUE <- FALSE
            if(sum(sort(snps) - snps)==0){
                SORT <- TRUE
            }
            if(length(unique(snps))-length(snps)==0){
                UNIQUE <- TRUE
            }
            return(data.frame(SORT=SORT, UNIQUE=UNIQUE))
            })
    if(FALSE %in% sort_unique$SORT){
        stop("snps must be sorted for each unqiue tetrad and chromosome combo")
    }
    if(FALSE %in% sort_unique$UNIQUE){
        stop("snps must be uniquely identified for each unqiue tetrad and chromosome combo")
    }   
}

# If direction is left (right), return the snp data and extent of the unique region 
#    immediately to the left (right) of the focal region.
get_flanking <- function(res_range, direction, df, threshold_size){
    if(direction=="left"){
        tmp <- df[which(df$snp==res_range[1])-1,]
    }
    if(direction=="right"){
        tmp <- df[which(df$snp==res_range[2])+1,]
    }
        flank_range <- range(CO_summary$snp[which(CO_summary$type==tmp$type)])
        flank_extent <- 1+(flank_range[2] - flank_range[1])
    if(tmp$GCbias==2){
        if(flank_extent<threshold_size){
            type2 <- "GC"  
        } else {
            type2 <- "nonGC"
        }
    } else {
        type2 <- "GC"
    }



    return(cbind(tmp, flank_extent, type2))
}

# A function for converting recombine object to tetrad.states
#' @export recombine_to_tetrad_states
recombine_to_tetrad_states <- function(tetrad_data){
    if(!inherits(tetrad_data, "recombine")){
        stop(paste("Object 'tetrad_data' needs to be of class 'recombine'.", sep=""))
    }
    n <- length(tetrad_data$parents$snps)
    out <- data.frame(Tetrad=rep(1, n), Chr=rep(1, n),
        Snp=tetrad_data$parents$snps, one=tetrad_data$chromatids.recombined$p1_1, 
        two=tetrad_data$chromatids.recombined$p1_2, three=tetrad_data$chromatids.recombined$p2_1, 
        four=tetrad_data$chromatids.recombined$p2_2)
    class(out) <- c("data.frame", "tetrad.states")
    return(out)

}

# Internal function that converts forward.backward class to tetrad.states class
fb_2_tetrad_states <- function(data){
    res <- FALSE
    if(inherits(data, "forward.backward")){
        res <- TRUE
    }
    if(!res){
        stop("Object data needs to be of class forward.backward.")
    }
        out <- data.frame(snp=data[[1]]$snp.locations, one=data[[1]]$states_inferred,
            two=data[[2]]$states_inferred,three=data[[3]]$states_inferred,four=data[[4]]$states_inferred)
        class(out) <- c("data.frame", "tetrad.states")
        return(out)        
    }

# Simple function to check if there is a tie in a max arguement
# Returns 1 if true and 0 if false (no tie)
check_tie <- function(x){
    # Do they equal each other? 
    if(length(which(x==max(x)))>1){
        return(1)
    }
    # Is there a unique maximum value? 
    if(length(which(x==max(x)))==1){
        return(0)
    }   
    # See note in forward-backward.R (section #6) for how this can occur. 
    if(NaN %in% x){
        return(1)
    } 
}

# Estimate crossover and gene conversion events
unique_regions <- function(data, threshold_size){

    # Initialize values
    sums <- apply(data[,4:7],1,sum)
    text <- apply(data[,4:7],1,function(x){paste(x[1],x[2],x[3],x[4],sep="_")})
    type <- numeric(length(sums))
    GC_count <- 0
    CO_count <- 0

    # Unique + values are sep recombination regions; unique - values are unique GC events
    for(a in 1:length(text)){
        # First position:      
        if(a==1){
            # Is this a GC?
            if(sums[a]!=2){
                GC_count <- GC_count - 1
                type[a] <- GC_count
            }  else {
                CO_count <- CO_count + 1
                type[a] <- CO_count
            } 
        }
        # All other positions:
        if(a>1){
            # If this is a GC:
            if(sums[a]!=2){
                if(text[a]==text[a-1]){
                    type[a] <- GC_count
                }
                if(text[a]!=text[a-1]){
                    GC_count <- GC_count - 1                        
                    type[a] <- GC_count
                }
            }
            # If this is a 2:2 ratio:
            if(sums[a]==2){
                if(text[a]!=text[a-1]){
                    CO_count <- CO_count + 1
                    type[a] <- CO_count
                }
                if(text[a]==text[a-1]){
                    type[a] <- CO_count
                }
            }
        }
    }

    out <- list(tetrad_states=data, type=type)
    class(out) <- c("list", "unique.regions")
    return(out)
}


check_GC_bias <- function(tetrad_states){
    if(!inherits(tetrad_states, "data.frame")){
        stop(paste("Object 'tetrad_states' needs to be of class 'data.frame'.", sep=""))
    }
    sums <- apply(tetrad_states[,4:7],1,sum)
    return(sums)
}


# This internal function returns the actual segregation pattern:
get_text <- function(i){
    return(paste(i[,'one'],i[,'two'],i[,'three'], i[,'four'], sep="_"))
}

# Infer additional tracts.
# For example, 'COnoGC' is a cross-over region where no gene conversion was detected. Because there are 0
# snps detecting gene conversion, the 'snp' location of the tracts is in between the crossover region.  
infer_COnoGC_tracts <- function(inferred_tracts){
    if(!inherits(inferred_tracts, "inferred.tracts")){
        stop(paste("Object 'inferred_tracts' needs to be of class 'inferred.tracts'.", sep=""))
    }
    tetrad <- inferred_tracts$tetrad[1]
    chr <- inferred_tracts$chr[1]

    # Call COnoGCs if present:
    COnoGCs <- lapply(1:max(1,(dim(inferred_tracts)[1]-1)), function(i, ...){
            # Check to see if the i-th row is a 2:2 tract
            if(inferred_tracts[i,'type']=="2_2" & dim(inferred_tracts)[1] >= (i+1)){
                # If so, check to see if the i+1 row is also a 2:2 ratio
                if(inferred_tracts[(i+1),'type']=="2_2"){
                    # If so, return a COnoGC entry
                    res_range <- c(inferred_tracts[i, 'end_snp'],inferred_tracts[(i+1), 'start_snp'])
                    pt <- mean(res_range)
                    return(data.frame(region=0, tetrad=tetrad, chr=chr,type="COnoGC", start_snp=pt, 
                        end_snp=pt, extent=res_range[2]-res_range[1], bias=NA))
                }
            }

        })
    COnoGCs <- do.call(rbind, COnoGCs)

    # Combine datasets and sort:
    out <- rbind(inferred_tracts, COnoGCs)
    out2 <- out[with(out, order(start_snp)),]

    # Return processed output:
    return(out2)

}


# merge_tracts <- function(inferred_tracts){
#     if(!inherits(inferred_tracts, "inferred.tracts")){
#         stop(paste("Object 'inferred_tracts' needs to be of class 'inferred.tracts'.", sep=""))
#     }
#     tetrad <- inferred_tracts$tetrad[1]
#     chr <- inferred_tracts$chr[1]



# }




# Checks if object i is of class tetrad.states or a list of 4 elements, each of 
#    class forward.backward
check_class <- function(i){
    res <- FALSE
    if(inherits(i, "tetrad.states") | inherits(i, "forward.backward")){
        res <- TRUE
    }
    return(res)
}

#' @export fwd_back_to_tetrad_states
fwd_back_to_tetrad_states <- function(fb_data){
    trim <- fb_data[,c('Tetrad', 'Spore', 'Chr', 'Snp', 'states_inferred')]
    trim1 <- reshape(trim, timevar="Spore", idvar=c("Tetrad", "Chr", "Snp"), direction="wide")
    colnames(trim1) <- c("Tetrad", "Chr", "Snp", "one", "two", "three", "four")
    class(trim1) <- c("data.frame", "tetrad.states")
    return(trim1)
}
