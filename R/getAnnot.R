#' Retrieve annotation information 
#' 
#' Uses the \code{annotatr} package to retrieve annotation information (
#' CpG category and gene coding sequences) for the \code{annoTrack} argument
#' of \code{\link{plotDMRs}}. Allows for 5 
#' re-tries if download fails (to allow for a spotty internet connection).
#' 
#' @details Note that this package needs to attach the \code{annotatr} package, 
#' and will
#' return NULL if this cannot be done. You can still use the 
#' \code{\link{plotDMRs}} function without this optional annotation step, 
#' just by leaving the \code{annoTrack} argument as NULL.
#' 
#' @param genomeName a character object that indicates which organism is 
#' under study. Use the function \code{builtin_genomes()} to see
#' a character vector of available genome names to choose from (see 
#' \code{annotatr} documentation for more details).
#' 
#' @return a \code{SimpleGRangesList} object with two elements returned
#' by \code{\link{getAnnot}}. The first
#' contains CpG category information in the first element (optional)
#' coding gene sequence information in the second element (optional).
#' At least one of these elements needs to be non-null in order for 
#' any annotation to be plotted, but it is not necessary to contain
#' both.
#' 
#' @export
#' 
#' @import annotatr
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom rtracklayer liftOver
#' @importFrom GenomeInfoDb genome
#' 
#' @examples
#' 
#' # get annotation information for hg19
#' annoTrack <- getAnnot('hg19')
#' 
#' 
getAnnot <- function(genomeName) {
    requireNamespace("annotatr")
    liftTo <- NULL
    if(genomeName == 'hg18'){
        message("Genome ", genomeName, " will be built by lifting over hg19 annotations")
        liftTo <- 'hg18'
        genomeName <- 'hg19'
    } else if (!genomeName %in% c(annotatr::builtin_genomes(), "Tthymallus", "ThyArc1.0")) {
        message("Genome ", genomeName, " is not supported by annotatr or custom Grayling scripts")
        return(NULL)
    }
    
    if (is.null(genomeName)) {
        return(NULL)
    } else {
        if (genomeName %in% c("Tthymallus", "ThyArc1.0")) {
            
            message(paste("Retrieving local annotations for", genomeName))
            cpg <- try(getCpGs(genome = genomeName), silent = FALSE)
            genes <- try(getGenes(genome = genomeName), silent = FALSE)      
            fail1 <- ifelse(is(cpg, "try-error"), 1, 0)
            fail2 <- ifelse(is(genes, "try-error"), 1, 0)
            
        } else {
            # 3. Standard annotatr Workflow (UCSC Genomes)
            annot_CpG <- paste0(c(genomeName, "_cpgs"), collapse = "")
            annot_genes <- paste0(c(genomeName, "_genes_cds"), collapse = "")
            
        }
        
        # 4. Processing and cleanup
        if (fail1 == 0 && fail2 == 0) {
            
            # Custom substrate string cleanup for T. arcticus
            if(genomeName %in% c("Dpulex", "Tthymallus", "ThyArc1.0")){
                # but if it does, keep it:
                # cpg$type <- substr(cpg$type, 1, nchar(cpg$type)) 
            } else {
                cpg$type <- substr(cpg$type, 10, nchar(cpg$type))
            }
            
            annot <- GRangesList(CpGs = cpg, Exons = genes, compress=FALSE)
            return(annot)
        } else {
            message("Annotation retrieval failed.")
            return(NULL)
        }
    }
}
