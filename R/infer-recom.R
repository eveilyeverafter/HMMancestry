#' @title Identify crossover, non-crossover, and gene conversion tracts
#' 
#' @description Infers different types of tracts (crossover, non-crossover, and gene conversion) along a 
#' chromosome
#' 
#' @param \code{data} A data.frame with 7 columns:
#'  \enumerate{
#'      \item{c("Tetrad", "Chr", "Snp", "one", "two", "three", "four")}
#'      \item{\code{Tetrad} specifying the tetrad ID}
#'      \item{\code{Chr} giving the chromosome name}
#'      \item{\code{Snp} a vector of snp locations (in bps)}
#'      \item{\code{one} the inferred states for spore 1}
#'      \item{\code{two} the inferred states for spore 2}
#'      \item{\code{three} the inferred states for spore 3}
#'      \item{\code{four} the inferred states for spore 4}
#' }
#' 
#' @param \code{threshold_size} The size (in bps) of the threshold (see details)
#' 
#' @details Uses the 3 step classification scheme described in Hether et al. (in review)
#' to identify the location of these specific CO, NCO, and telomeric gene conversion tracts.
#' Specifically, \code{infer_tracts} attempts to identify the following tracts:
#'      \itemize{
#'          \item{2_2}{ a tract with 2:2 segregation}
#'          \item{COyesGC}{ a region of gene conversion that was associated with a crossover event}
#'          \item{COnoGC}{ a region between a crossover event that did not have a detectable gene conversion}
#'          \item{NCO}{ a non-crossover; a gene conversion tract without crossing}
#'          \item{GC_tel}{ a gene conversion tract located at the chromosome end}
#'      }  
#'
#' @return A data.frame containing the following columns:
#' \enumerate{
#'      \item{type}{ The type of inferred tract}
#'      \item{start_snp}{ The starting snp position in base pairs}
#'      \item{end_snp}{ The ending snp position in base pairs}
#'      \item{extent}{ The size of the tract. For COnoGC, this extent is the s
#'            spanning region between flanking CO events.}
#'      } 
#' 
#' @author Tyler D. Hether
#' 
#' @export infer_tracts
#'
#' @references Hether, T.D., C. G. Wiench1, and P.A. Hohenlohe (in review). 2015. Novel molecular and analytical tools 
#' for efficient estimation of rates of meiotic crossover, non-crossover and gene conversion
#' 
#' @examples
#' set.seed(1234567)        # For reproducibility
#' n_tetrads <- 3           # number of tetrads
#' l <- 1000                # number of snps to simulate
#' c <- 3e-05               # recombination rate between snps (Morgan/bp)
#' snps <- c(1:l)*250       # snps are evenly spaced 250 bp apart
#' p_a <- 0.95              # assignment probability
#' coverage <- 1            # mean coverage
#' # simulate tetrads
#' tetrad <- sim_tetrad(n.tetrads=n_tetrads, scale=c, snps=snps, 
#'  p.assign=p_a, mu.rate=1e-03, f.cross=0.6, f.convert=0.8, 
#'  length.conversion=2e3, coverage=coverage)
#' #' # Example 1 -- infer tracts directly from simulated data
#' inf_tracts_sim <- infer_tracts(tetrad)
#' inf_tracts_sim
#' #' # Example 2 -- infer tracts from inferred data
#' inf_states <- ddply(tetrad, .(Tetrad, Spore, Chr), 
#'     function(x){
#'     return(fb_haploid(snp_locations=x$Snp, p0=x$p0, 
#'     p1=x$p1, p_assign=p_a, scale=c))})
#' inf_tracts_inf_states <- infer_tracts(inf_states)
#' inf_tracts_inf_states

infer_tracts <- function(dat, threshold_size=2.5e3){
    dat <- as.data.frame(dat)
    if(dim(dat)[2]!=7 & dim(dat)[2]!=18){
        stop(paste("object dat needs to come from either sim_tetrad or fb_haploid"))
    }
    if(dim(dat)[2]==18){
        dat <- dat[,c(1,2,3,4,5,6,17)]
    }
    colnames(dat) <- c("Tetrad", "Spore", "Chr",  "Snp", "p0", "p1", "states")

    
    out_all <- plyr::ddply(as.data.frame(dat), .(Tetrad, Chr), function(data){

        # Reshape from long to wide
        reshaped_data <- fb_to_tetrad_states(data)
        # Check and sort the data
        checked_data <- prep_infer_tracts_data(reshaped_data) 

        # identify regions in each tetrad and each chromosome
        CO_summary <- summarize_regions(checked_data)
        # Now do the algorithm (both parts)
        # What are the unique regions?
        regions <- unique(CO_summary$type)

        outlist <- lapply(regions, function(res, ...){
        
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
            return(data.frame(res, tetrad=CO_summary$Tetrad[1], chr=CO_summary$Chr[1], type, 
                res_range[1],res_range[2], extent,BIAS, pattern=as.character(CO_summary[CO_summary$type==unique(CO_summary$type)[pos],'text'][1])))
        })
        

        out2 <- do.call(rbind, outlist)




        # Return output of inferred tracts along the chromosome 
        # out <- out[,-1] # Housekeeping
        colnames(out2) <- c("region", "tetrad", "chr", "type", "start_snp", "end_snp", "extent", "bias", "pattern")
        # class(out2) <- c("data.frame", "inferred.tracts")
        class(out2) <- c("data.frame")


        # Reclassifying temporary tracts:


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

        # This identifies the approx location of COs without a detected gene conversion tract
        # print("Screening for COs without GCs...")
        out4 <- infer_COnoGC_tracts(out3)

        # Housekeeping:
        out4 <- out4[,-c(1,2,3)]
        return(out4)
    })
   
    # Housekeeping: 
    # OUT <- OUT[,-c(3,4,5)]
    return(out_all)
}

#' @title Identify crossover points from state sequences and snp locations
#' 
#' @description This is a simple function that takes a vector of parental states
#' along a chromosome and identifies where states have changed.  
#' 
#' @param \code{snps_genotypes_df} a data.frame with two columns:
#'   \itemize{
#'      \item{snps}{ a vector of ordered snp locations (in bps)}
#'      \item{states}{ a vector of corresponding inferred states, either haploid
#'       (2 states = 0 or 1) or diploid with three states possible (0,1,2)}
#'    }
#'
#' @return a data.frame containing the midpoint between focal snp \eqn{i} and snp 
#' \eqn{i-1} and whether their states were the same (0) or different (1). 
#' 
#' @seealso \code{\link{fb_haploid}}, \code{\link{fb_diploid}}
#'
#' @author Tyler D. Hether 
#' 
#' @export id_recombination_events
#' 
#' @examples
#' df <- data.frame(snps=100*(1:10), 
#'    states=c(rep(0,4), rep(1,6)))
#' id_recombination_events(df)

id_recombination_events <- function(snps_genotypes_df){ 

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
        flank_range <- range(df$snp[which(df$type==tmp$type)])
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

# A function for converting a haploid fb dataset to a tetrad_states (used in infer_tracts)
#' @export fb_to_tetrad_states
fb_to_tetrad_states <- function(fb_data){
    # For input, take the output of sim_tetrad, fb_haploid
    if(dim(fb_data)[2]==7){
        TMP <- fb_data[,c(1,2,3,4,7)]
       
    } else {
        TMP <- fb_data[,c('Tetrad', 'Spore', 'Chr', 'Snp', 'states_inferred')]
    }
    colnames(TMP) <- c("Tetrad", "Spore", "Chr", "Snp", "states")
    
    w <- reshape(TMP, 
         timevar = "Spore",
         idvar = c("Tetrad", "Chr", "Snp"),
         direction = "wide")

    if(!all(dim(w[complete.cases(w),])==dim(w))){
        warning(paste("Some spores had missing data at 1 or more snps in Tetrad ", TMP[1,1], ", Chr ", TMP[1,2], ". Those snps will be removed.", sep=""))
        w <- w[complete.cases(w),]
    }   
    colnames(w) <- c("Tetrad", "Chr", "Snp", "one", "two", "three", "four")
    return(w)
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
    # if(!inherits(inferred_tracts, "inferred.tracts")){
    #     stop(paste("Object 'inferred_tracts' needs to be of class 'inferred.tracts'.", sep=""))
    # }
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
                        end_snp=pt, extent=res_range[2]-res_range[1]))#, bias=NA, pattern="9_9_9_9"))
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

# Internal functiont that checks the data format and sort it if necessary.
prep_infer_tracts_data <- function(data){
    if(!("Tetrad" %in% colnames(data))){
        stop("No Tetrad column detected")
         # inferred_states_data$Tetrad <- as.numeric(unlist(lapply(strsplit(inferred_states_data$Ind, split="_"), "[", 1L)))
    } 
    # data is a data.frame with the following columns
    colnames(data) <- c("Tetrad", "Chr", "Snp", "one", "two", "three", "four")
    # Order the data
    data2 <- data[with(data, order(Tetrad, Chr, Snp)),]
    # filter out any NaNs
    data2 <- data[complete.cases(data), ]
    return(data2)
}
# Internal function to call unique regions for each tetrad and each chromosome
summarize_regions <- function(dat){
    out <- ddply(dat, .(Tetrad, Chr), function(tetrad_states){
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
        to_return <- data.frame(snp=regions$tetrad_states[,'Snp'], type=regions$type, GCbias=biases, text=text)
        return(to_return)
    })
    return(out)
}

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




