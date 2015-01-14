#' @title Identify crossover, non-crossover, and gene conversion trakcs
#' 
#' @description Infers different types of tracks (crossover, non-crossover, and gene conversion) along a chromosome given genotyped yeast tetrads or simulated data
#' 
#' @param tetrad_states a \code{states.matrix} object inherited from either \code{XXX} or \code{YYY} 
#' 
#' @param tetrad (Optional) A numeric or character specifying the tetrad ID
#' 
#' @param chr (Optional) A numeric or character specifying the chromosome ID
#' 
#' @details (to do)
#' 
#' @return A data.frame containing the following columns:
#' \describe{
#'     \item{tetrad}{The tetrad ID (default == 1)}
#'     \item{chr}{The chromosome ID (default == "I")}
#'     \item{type}{The type of inferred track:
#'         \describe{
#'             \item{2_2}{a track with a 2:2 ratio}
#'             \item{GC_tel}{a gene conversion track located at a chromosomal end}
#'             \item{GC_internal}{a gene conversion track with at least one flanking track also a gene conversion track}
#'             \item{COnoGC}{indicates where a crossover event occurred but without a detectable gene conversion track}
#'             \item{COyesGC}{a gene conversion track associated with a crossover event}
#'             \item{NCO}{a gene conversion track without an associated crossover event}
#'             }}
#'     \item{start_snp}{numeric value of the first snp in a given track}
#'     \item{end_snp}{numeric value of the last snp in a given track}
#'     \item{extent}{numeric value of the range of snps (inclusive) in a given track}
#'     \item{bias}{integer specifying the gene conversion bias where:
#'         \describe{
#'         \item{0}{a 4:0 segregation (all of parent 0 type)
#'         \item{1}{a 3:1 segregation
#'         \item{2}{a 2:2 no bias found
#'         \item{3}{a 1:3 segregation
#'         \item{4}{a 0:4 segregation (all of parent 1 type)
#'         \item{NA}{(COnoGC only) bias information is not applicable}    
#'         }}
#' }
#' 
#' @seealso \code{recombine_index}, \code{recombine}, \code{recombine_to_tetrad_states}
#' 
#' @author Tyler D. Hether
#' 
#' @export infer_tracks
#' 
#' @examples
# Simulated example 1
set.seed(1234567) # For reproducability
l <- 1000 # number of loci to simulate
rec <- 0.01 # recombination rate between each snp
r <- recombine_index(rep(rec, l-1)) # recombination rate between each snp (vector form)
p_a <- .999 # probability of correct sequencing assignment (1-sequence error rate)
p <- make_parents(l) # make the parent
recomb_sim <- recombine(parents=p, r.index=r, mu.rate=0) # recombine parents
states <- recombine_to_tetrad_states(tetrad_data=recomb_sim) # convert to tetrad.states object

# add a couple of gene conversion events
states[1:10,'three'] <- 1 # GC_tel
states[1:4,'two'] <- 1 # GC_internal
states[847:855,'three'] <- 1 # COyesGC
states[944:971,'one'] <- 0 # NCO

infer_tracks(tetrad_states=states, tetrad=1, chr="I")

# Simulated example 2, using the fb algorithm on 1x sequencing coverage
set.seed(1234567) # For reproducability
l <- 1000 # number of loci to simulate
rec <- 0.01 # recombination rate between each snp
r <- recombine_index(rep(rec, l-1)) # recombination rate between each snp (vector form)
p_a <- .999 # probability of correct sequencing assignment (1-sequence error rate)
p <- make_parents(l) # make the parent
recomb_sim <- recombine(parents=p, r.index=r, mu.rate=0) # recombine parents
sim_reads <- simulate_coverage(a=recomb_sim, p_assign=p_a, coverage=1) # simulate sequencing coverage


# Use the forward-backward algorithm to get the posterior probability of parent '0' ancestry
fbres <- lapply(c(1:4), function(x,...){
    out <- estimate_anc_fwd_back(snp_dat=sim_reads, spore_number=x, 
    chr_name="I", p_assign=p_a, p_trans=rec)
    return(out)
    })




# # Convert to matrix and add snp names:
# tetrad_states <- cbind(fbres[[1]]$snp_locations, matrix(unlist(state_res), ncol=4, byrow=FALSE))
# colnames(tetrad_states) <- c("snp", "one", "two", "three", "four")
# # Remove states with ambiguous calls (i.e., PP=0.5): 
# states <- tetrad_states[complete.cases(tetrad_states),]
# ~~~~~~~~~~~~~~~~~~~






# Main function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


infer_tracks <- function(tetrad_states=NULL, forward_backward_data=NULL, tetrad=1, chr="I"){
    # Two datasets are allowed to pass through this function.
    # 1) of class "tetrad.states" generally inherited from the recombine_to_tetrad_states function OR
    # 2) a list of four elements, each of class forward.backward that is a result of the estimate_anc_fwd_back function.
FIX THIS BLOCK   
    if(tetrad_states!=NULL){
        if(!inherits(tetrad_states, "tetrad.states")){
            stop(paste("Object 'tetrad_states' needs to be of class 'data.frame'.", sep=""))
        }
    }

     if(forward_backward_data!=NULL){
        if(!inherits(tetrad_states, "forward.backward")){
            stop(paste("Object 'forward_backward_data' needs to be of class 'forward.backward'.", sep=""))
        }
    }
FIX THAT BLOCK   

    # Get the unique candidate 'regions':
    regions <- unique_regions(tetrad_states)

    # Calculate the biases, if present:
    biases <- check_GC_bias(tetrad_states)

    # Get the actual segregation pattern:
    text <- get_text(tetrad_states)

    # Combine datasets: 
    CO_summary <- data.frame(snp=regions$tetrad_states[,'snp'], type=regions$type, GCbias=biases, text=text)

    # Go through the chromosome and call each region as a specific 'track':
    out <-   do.call(rbind, lapply(unique(CO_summary$type), function(res, ...){
        type="NULL"
        BIAS="NULL"
        # print(res)
        res_range <- range(CO_summary$snp[which(CO_summary$type==res)])
        extent <- 1+(res_range[2] - res_range[1])
        
        # Call a 2:2 region if present:
        if(res>0){
            type <- "2_2"
            BIAS <- 2
        }
        
        # If the unique region is not a 2:2 region, then do the following:  
        if(res<0){
            pos <- which(unique(CO_summary$type)==res)
            
            # If it's a non-2:2 region and it's located on the chromosomal end, it's not internal:
            if(res==unique(CO_summary$type)[1] | res==rev(unique(CO_summary$type))[1]){
                internal <- FALSE
            } else {
                internal <- TRUE
            }

            # If this non 2:2 region is located on the end of the chromosome, call it a GC_tel:
            if(!internal){                        
                type <- "GC_tel"
            }
            # Is this non 2:2 region located internally?
            if(internal){
                # First determine if there are 2:2 regions immediately flanking this region.
                if(CO_summary$GCbias[res_range[1]-1]==2 & CO_summary$GCbias[res_range[2]+1]==2){
                        # If so, call either a NCO if flanking regions are identical or a GCwCO if 
                        # the flanking regions are not identical.
                        if(as.character(CO_summary$text[res_range[1]-1]) == as.character(CO_summary$text[res_range[2]+1])){
                            type <- "NCO"
                        } else {
                            type <- "COyesGC"
                        }

                }

                # Second, if this internal region isn't flanked by two 2:2 regions, then
                # call it a GC_internal to be called later. 
                if(CO_summary$GCbias[res_range[1]-1]!=2 | CO_summary$GCbias[res_range[2]+1]!=2){
                    type <- "GC_internal"
                }

            }    

            # Last, return the ancestry bias
            BIAS <- CO_summary[CO_summary$type==unique(CO_summary$type)[pos],'GCbias'][1]
        }

        # Store results as a data.frame
        return(data.frame(res, tetrad, chr, type, res_range[1],res_range[2], extent,BIAS))
    }))
    
    # Return output of inferred tracks along the chromosome 
    out <- out[,-1] # Housekeeping
    colnames(out) <- c("tetrad", "chr", "type", "start_snp", "end_snp", "extent", "bias")
    class(out) <- c("data.frame", "inferred.tracks")

    # Call additional tracks:
    out2 <- infer_COnoGC_tracks(out)
    return(out2)
}




# Minor function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# A function for converting recombine object to tetrad.states
recombine_to_tetrad_states <- function(tetrad_data){
    if(!inherits(tetrad_data, "recombine")){
        stop(paste("Object 'tetrad_data' needs to be of class 'recombine'.", sep=""))
    }

    out <- data.frame(snp=tetrad_data$parents$snps, one=tetrad_data$chromotids.recombined$p1_1, two=tetrad_data$chromotids.recombined$p1_2, 
    three=tetrad_data$chromotids.recombined$p2_1, four=tetrad_data$chromotids.recombined$p2_2)
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
unique_regions <- function(tetrad_states){
    if(!inherits(tetrad_states, "data.frame")){
        stop(paste("Object 'tetrad_states' needs to be of class 'data.frame'.", sep=""))
    }

    # Initialize values
    sums <- apply(tetrad_states[,2:5],1,sum)
    text <- apply(tetrad_states[,2:5],1,function(x){paste(x[1],x[2],x[3],x[4],sep="_")})
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

    out <- list(tetrad_states=tetrad_states, type=type)
    class(out) <- c("list", "unique.regions")
    return(out)
}


check_GC_bias <- function(tetrad_states){
    if(!inherits(tetrad_states, "data.frame")){
        stop(paste("Object 'tetrad_states' needs to be of class 'data.frame'.", sep=""))
    }
    sums <- apply(tetrad_states[,2:5],1,sum)
    return(sums)
}


# This internal function returns the actual segregation pattern:
get_text <- function(tetrad_states){
    if(!inherits(tetrad_states, "data.frame")){
        stop(paste("Object 'tetrad_states' needs to be of class 'data.frame'.", sep=""))
    }

    return(paste(tetrad_states[,'one'],tetrad_states[,'two'],tetrad_states[,'three'], tetrad_states[,'four'], sep="_"))
}

# Infer additional tracks.
# For example, 'COnoGC' is a cross-over region where no gene conversion was detected. Because there are 0
# snps detecting gene conversion, the 'snp' location of the tracks is in between the crossover region.  
infer_COnoGC_tracks <- function(inferred_tracks){
    if(!inherits(inferred_tracks, "inferred.tracks")){
        stop(paste("Object 'inferred_tracks' needs to be of class 'inferred.tracks'.", sep=""))
    }
    tetrad <- inferred_tracks$tetrad[1]
    chr <- inferred_tracks$chr[1]

    # Call COnoGCs if present:
    COnoGCs <- lapply(1:(dim(inferred_tracks)[1]-1), function(i, ...){
            # Check to see if the i-th row is a 2:2 track
            if(inferred_tracks[i,'type']=="2_2"){
                # If so, check to see if the i+1 row is also a 2:2 ratio
                if(inferred_tracks[(i+1),'type']=="2_2"){
                    # If so, return a COnoGC entry
                    res_range <- c(inferred_tracks[i, 'end_snp'],inferred_tracks[(i+1), 'start_snp'])
                    pt <- mean(res_range)
                    return(data.frame(tetrad=tetrad, chr=chr,type="COnoGC", start_snp=pt, end_snp=pt, extent=res_range[2]-res_range[1], bias=NA))
                }
            }

        })
    COnoGCs <- do.call(rbind, COnoGCs)

    # Combine datasets and sort:
    out <- rbind(inferred_tracks, COnoGCs)
    out2 <- out[with(out, order(start_snp)),]

    # Return processed output:
    return(out2)

}



















