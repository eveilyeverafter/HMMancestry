#' Function for building data
#'
#' @param snp_spore_1 A list of four elements, each containing a named vector of SNPs, indicating number of reads for each snp
make_snp_data <- function(parent_1, parent_2, snp_id){

    dat <- check_snp_calls(parent_1, parent_2, snp_id)
    out <- build_snp_table(dat)
    class(out) <- c("data.frame", "snp.recom")
    out  
}

check_snp_calls <- function(parent_1, parent_2, snp_id){

    if (!inherit(parent_1, "list") | !inherit(parent_2, "list"))
        stop("Both 'snp_spore_1' and 'snp_spore_2' must be lists")

    if (length(parent_1) != 4)
        stop("parent_1 must have four elements")

    if (length(parent_2) != 4)
        stop("parent_2 must have four elements")

    if (names(parent_1) != names(parent_2))
        stop("snp names need to match")

    if (!all(is.integer(parent_1)) | !all(is.integer(parent_2)))
        stop("read counts must be integers for all SNPs")

    if (!all(names(parent_1)) %in% names(snp_names))
        stop("names of snp_id need to match those of snps in parents")
    
    snp_id <- snp_id[snp_names]

    list(parent_1=parent_1, parent_2=parent_2, snp_id=snp_id)

}

build_snp_table <- function(dat){

    p1 <- dat$parent_1
    p2 <- dat$parent_2
    id <- dat$snp_id

    ## build whatever table you want
    
}
    

    
