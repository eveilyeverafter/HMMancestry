#' @title Infer crossover, non-crossover, and gene conversion events.
#'
#' @description Infers the location and type of each recombination 'track'. 


set.seed(1234567)
l <- 1000
rec <- 0.01
p_a <- .999
p <- make_parents(l)
r <- recombine_index(rep(rec, l-1))
a <- recombine(parents=p, r.index=r, mu.rate=0)
sim_reads <- simulate_coverage(a=a, p_assign=.999, coverage=10)

# Simulate sum data
fbres <- lapply(c(1:4), function(x,...){
    out <- estimate_anc_fwd_back(snp_dat=sim_reads, spore_number=x, 
    chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)
    return(out)
    })

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

# From the posterior probs, call states 0 or 1. 
# In the rare case of a tie (PP=0.5 for both), flag that state as -1.
state_res <- lapply(1:4, function(x){
      apply(fbres[[x]]$posterior, 1, function(i){
        if(NaN %in% i){
            warning("Something might not be right. NaNs are present in data.")
        }
        if(check_tie(i)==0){
            return(which(i==max(i))-1)
        }
        if(check_tie(i)==1){
            return(NA)
        }
      })
    })

# Convert to matrix and add snp names:
states_mat <- cbind(fbres[[1]]$snp_locations, matrix(unlist(state_res), ncol=4, byrow=FALSE))
colnames(states_mat) <- c("snp", "one", "two", "three", "four")
# Remove states with ambiguous calls (i.e., PP=0.5): 
states <- states_mat[complete.cases(states_mat),]

# ^ ^ ^ So now we have a snp-by-5 matrix specifying the location of snps (first column)
# and the ancestry assignment for each of the four spores (columns 2:5). 
# The next step is to summarize the data into unique 'events' -- regions of variable length
# with the the same recombination tract (i.e., 2:2 tract, gene conversion tract, etc.)

 
# The current simulation code ignores gene conversions, GCs. 
# In cases where there is low sequencing coverage, the FB may take time to switch from one
# parent to the other and this effect looks like a gene conversion event. 
# 
# Simulate telomeric gene conversion:
states[1:10,3] <- 1 # GC_tel
states[1:4,2] <- 1 # GC_tel

class(states) <- c("matrix", "states.matrix")

# Estimate crossover and gene conversion events
unique_regions <- function(states_mat){
    if(!inherits(states_mat, "states.matrix")){
        stop(paste("Object 'states_mat' needs to be of class 'states.matrix'.", sep=""))
    }

    # Initialize values
    sums <- apply(states_mat[,2:5],1,sum)
    text <- apply(states_mat[,2:5],1,function(x){paste(x[1],x[2],x[3],x[4],sep="_")})
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

    out <- list(states_mat=states_mat, type=type)
    class(out) <- c("list", "unique.regions")
    return(out)
}

# unique_regions(states_mat=states)



# For all regions lacking a 2:2 segregation, return the ancestral bias (0, or 1)
check_GC_bias <- function(states_mat){
    if(!inherits(states_mat, "states.matrix")){
        stop(paste("Object 'states_mat' needs to be of class 'states.matrix'.", sep=""))
    }
    sums <- apply(states_mat[,2:5],1,sum)
    # 0 = bias towards parent 0; 1 = bias towards parent 1; -1 = no bias (2:2).
    out <- data.frame(snps=states_mat[,1], bias=sapply(sums, function(x){
            if(x<2){
                return(0)
            }
            if(x==2){
                return(-1)
            }
            if(x>2){
                return(1)
            }
            }))
    # out <- all[complete.cases(all),]
    return(out)
}

# check_GC_bias(states_mat=states)

# This internal function returns the actual segregation pattern:
get_text <- function(states_mat){
    if(!inherits(states_mat, "states.matrix")){
        stop(paste("Object 'states_mat' needs to be of class 'states.matrix'.", sep=""))
    }

    return(paste(states_mat[,'one'],states_mat[,'two'],states_mat[,'three'], states_mat[,'four'], sep="_"))
}


# Infers the type of track along each position of a chromosome:
infer_tracks <- function(states_mat, chr="I"){
    if(!inherits(states_mat, "states.matrix")){
        stop(paste("Object 'states_mat' needs to be of class 'states.matrix'.", sep=""))
    }

    # Get the unique candidate 'regions':
    regions <- unique_regions(states_mat)

    # Determine if there is a bias (non-2:2 segregation):
    biases <- check_GC_bias(states_mat)

    # Get the actual segregation pattern:
    text <- get_text(states_mat)

    # Combine datasets: 
    CO_summary <- data.frame(snp=regions$states_mat[,'snp'], type=regions$type, GCbias=biases$bias, text=text)

    # Go through the chromosome and call each region as a specific 'track':
    # for(res in unique(CO_summary$type)){
    infer <-   sapply(unique(CO_summary$type), function(res, ...){
        type="NULL"
        print(res)
        res_range <- range(CO_summary$snp[which(CO_summary$type==res)])
        extent <- res_range[2] - res_range[1]
        
        # Call a 2:2 region if present:
        if(res>0){
            type="2_2"
            BIAS=-1
        }
         
        if(res<0){
            pos <- which(unique(CO_summary$type)==res)
            
            # If it's a non-2:2 region and it's located on the chromosomal end, it's not internal:
            if(pos==1 | pos==rev(unique(CO_summary$type))[1]){
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
                if(CO_summary$GCbias[res_range[1]+1]==0 & CO_summary$GCbias[res_range[2]+1]==0){
                        # If so, call either a NCO if flanking regions are identical or a GCwCO if 
                        # the flanking regions are not identical.
                        if(as.character(CO_summary$text[res_range[1]+1]) != as.character(CO_summary$text[res_range[2]+1])){
                            type <- "NCO"
                        } else {
                            type <- "GCwCO"
                        }

                }

                # Second, if this internal region isn't flanked by two 2:2 regions, then
                # call it a GC_unassigned to be called later. 
                if(CO_summary$GCbias[res_range[1]+1]!=0 | CO_summary$GCbias[res_range[2]+1]!=0){
                    type <- "GC_unassigned"
                }

            }    

            # Last, return the ancestry bias
            BIAS <- CO_summary[CO_summary$type==unique(CO_summary$type)[pos],'GCbias'][1]

        }

        system(paste('echo', paste(chr, type, res_range[1],res_range[2], extent,BIAS,sep=","), paste('>>./res_', chr, '.txt',sep=""), sep=" "))
        # return(data.frame(type=type, res_range1=res_range[1], res_range2=res_range[2], extent=extent, BIAS=BIAS))
        
    })

}







