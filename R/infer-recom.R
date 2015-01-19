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
#' Example 1: 1 tetrad
#' set.seed(1234567) # For reproducability
#' l <- 1000 # number of loci to simulate
#' rec <- 0.01 # recombination rate between each snp
#' r <- recombine_index(rep(rec, l-1)) # recombination rate between each snp (vector form)
#' p_a <- .999 # probability of correct sequencing assignment (1-sequence error rate)
#' p <- make_parents(l) # make the parent
#' recomb_sim <- recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert=1, length.conversion=10) # recombine parents
#' states <- recombine_to_tetrad_states(tetrad_data=recomb_sim) # convert to tetrad.states object
#' states
#' #
#' # Example 2: 100 simulated tetrads with a recombination hotspot
#' set.seed(1234567) # For reproducability
#' l <- 200 # loci
#' # simulate a recombination hotspot between the 99th and 100th snp
#' rec <- c(rep(0.001, 99), 0.4, rep(0.001, 99))
#' n.tetrads <- 500 # number of spores to simulate
#' res <- sim_tetrad(n.tetrads=n.tetrads, l=l, rec=rec, p.assign=0.999, 
#'    mu.rate=0, f.cross=0.8, f.convert=0, length.conversion=0, coverage=0.5)
#' # loop through the results list and estimate parental states for each of the n.tetrads
#' fblist <- lapply(1:length(res), function(Z){
#'         tetrad_est_fwd_back(snp.dat=res[[Z]], tetrad.id=Z, chr.name="I", 
#'     p.assign=0.999, p.trans=mean(rec))   
#'     })
#' # loop through the estimated states paths and infer recombination events:
#' inferlist <- lapply(1:length(fblist), function(Z){
#'     infer_tracts(data=fblist[[Z]], tetrad=Z)
#'     })
#' # combine results into a single data.frame
#' df <- do.call(rbind, inferlist)
#' hist(filter(df, type=="COyesGC" | type=="COnoGC" | type=="NCO")$start_snp,
#'  breaks=l, xlab="snp", main="start position of recombination point")

infer_tracts <- function(data, tetrad=1, chr="I"){
    # data can be an object of two classes.
    # 1) of class "tetrad.states" generally inherited from the recombine_to_tetrad_states function OR
    # 2) a list of four elements, each of class forward.backward that is a result of the estimate_anc_fwd_back function. 
    if(!check_class(data)){
            stop(paste("Object 'data' needs to be of class tetrad.states or forward.backward", sep=""))
    }

    # Convert from class foward.backward to tetrad.states, if needed:
    if(!inherits(data, "tetrad.states")){
        tetrad_states <- fb_2_tetrad_states(data)
        chr <- data[[1]]$chr.name
        tetrad <- data[[1]]$tetrad.id
    } else {
        tetrad_states <- data
    }

    # Get the unique candidate 'regions':
    regions <- unique_regions(tetrad_states)

    # Calculate the biases, if present:
    biases <- check_GC_bias(tetrad_states)

    # Get the actual segregation pattern:
    text <- get_text(tetrad_states)

    # Combine datasets: 
    CO_summary <- data.frame(snp=regions$tetrad_states[,'snp'], type=regions$type, GCbias=biases, text=text)

    # Go through the chromosome and call each region as a specific 'tract':
    out <-   do.call(rbind, lapply(unique(CO_summary$type), function(res, ...){
        # print(res)
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
                # First, identify the snps that flank left and right of res 
                left <- get_flanking(res_range=res_range, direction="left", df=CO_summary)
                right <- get_flanking(res_range=res_range, direction="right", df=CO_summary)                
                # Second determine if there are 2:2 regions immediately flanking this region.
                if(left$GCbias==2 & right$GCbias==2){
                        # If so, call either a NCO if flanking regions are identical or a GCwCO if 
                        # the flanking regions are not identical.
                        if(as.character(left$text) == as.character(right$text)){
                            type <- "NCO"
                        } else {
                            type <- "COyesGC"
                        }

                }

                # Third, if this internal region isn't flanked by two 2:2 regions, then
                # call it a GC_internal to be called later. 
                if(left$GCbias!=2 | right$GCbias!=2){
                    type <- "GC_internal"
                }

            }    

            # Last, return the ancestry bias
            BIAS <- CO_summary[CO_summary$type==unique(CO_summary$type)[pos],'GCbias'][1]
        }

        # Store results as a data.frame
        return(data.frame(res, tetrad, chr, type, res_range[1],res_range[2], extent,BIAS))
    }))
    
    # Return output of inferred tracts along the chromosome 
    out <- out[,-1] # Housekeeping
    colnames(out) <- c("tetrad", "chr", "type", "start_snp", "end_snp", "extent", "bias")
    class(out) <- c("data.frame", "inferred.tracts")

    # Call additional tracts:
    out2 <- infer_COnoGC_tracts(out)
    return(out2)
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
#' @export id_hotspots
#' 
#' @examples
#' # Example 1: simple
#' # A recombination occurred between snp 3 and 4 and between 8 and 9.
#' statepath <- c(0,0,0,1,1,1,1,1,0,0,0)
#' id_hotspots(state.vector=statepath)
#' which(id_hotspots(state.vector=statepath)==1)
#' #
#' # Example 2: complex
#' set.seed(1234567) # For reproducability
#' # simulate a recombination hotspot between the 99th and 100th snp
#' rec <- c(rep(0.001, 99), 0.1, rep(0.001, 99))
#' # simulate 500 spores en masse
#' n.spores <- 500 
#' spores <- sim_en_masse(n.spores=n.spores, l=200, rec=rec, 
#'  p.assign=.999, mu.rate=0.001, f.cross=0.5, 
#'     f.convert=0.5, length.conversion=10, coverage=1)
#' 
#' # Run the fb algorithm to estimate the parental states:
#' Allspores <- lapply(1:n.spores, function(Z){
#'         fbres2 <- est_fwd_back(single.snp.dat=spores[[Z]], 
#'             spore_number=Z, chr.name="I", p.assign=0.999, 
#'             p.trans=mean(rec))
#'         return(fbres2)
#'     })
#' # extract the inferred states
#' Allspores2 <- lapply(Allspores, function(i){
#'     return(as.numeric(i$states_inferred))})
#' 
#' # organize into a df to apply
#' df <- do.call(rbind,lapply(Allspores, function(i){
#'  return(as.numeric(i$states_inferred))}))
#' 
#' # plot where a recombination event was inferred
#' plot(apply(t(apply(df, 1, id_hotspots)), 2, sum), 
#'  type="l", xlab="snp", ylab="Number of recombination events")

id_hotspots <- function(state.vector){
    if(length(state.vector)<=1){
        stop("At least two snps needed to infer recombination points.")
    }
    out <- numeric(length(state.vector))
    
    for(i in 2:length(state.vector)){
        if(state.vector[i]!=0 & state.vector[i]!=1){
            stop("states need to be 0 or 1 only.")
        }
        if(state.vector[i]!=state.vector[(i-1)]){
            out[i] <- 1
        } 
    }
    return(out)
}

# Minor functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If direction is left (right), return the snp data of the unique region 
#    immediately to the left (right) of the focal region.
get_flanking <- function(res_range, direction, df){
    if(direction=="left"){
        return(df[which(df$snp==res_range[1]-1),])
    }
    if(direction=="right"){
        return(df[which(df$snp==res_range[2]+1),])
    }
}

# A function for converting recombine object to tetrad.states
recombine_to_tetrad_states <- function(tetrad_data){
    if(!inherits(tetrad_data, "recombine")){
        stop(paste("Object 'tetrad_data' needs to be of class 'recombine'.", sep=""))
    }

    out <- data.frame(snp=tetrad_data$parents$snps, one=tetrad_data$chromatids.recombined$p1_1, 
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
unique_regions <- function(data){
    if(!check_class(data)){
            stop(paste("Object 'data' needs to be of class tetrad.states or forward.backward", sep=""))
    }

    # Initialize values
    sums <- apply(data[,2:5],1,sum)
    text <- apply(data[,2:5],1,function(x){paste(x[1],x[2],x[3],x[4],sep="_")})
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
    sums <- apply(tetrad_states[,2:5],1,sum)
    return(sums)
}


# This internal function returns the actual segregation pattern:
get_text <- function(i){
    if(!check_class(i)){
            stop(paste("Object 'data' needs to be of class tetrad.states or forward.backward", sep=""))
    }
    if(inherits(i, "tetrad.states")){
        return(paste(i[,'one'],i[,'two'],i[,'three'], i[,'four'], sep="_"))
        } else {
        return(paste(i[[1]]$states_inferred, i[[2]]$states_inferred, 
            i[[3]]$states_inferred, i[[4]]$states_inferred,sep="_"))
    }
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
                    return(data.frame(tetrad=tetrad, chr=chr,type="COnoGC", start_snp=pt, 
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

# Checks if object i is of class tetrad.states or a list of 4 elements, each of 
#    class forward.backward
check_class <- function(i){
    res <- FALSE
    if(inherits(i, "tetrad.states") | inherits(i, "forward.backward")){
        res <- TRUE
    }
    return(res)
}

















