###
### Data preparation scripts for three-group CER
################################################################################

###
### Functions for PS estimation
################################################################################

###  Generalized PS Function
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @param data data_frame
##' @param formula formula for PS model
##' @param family default \code{multinomial(parallel = FALSE)}
##' @param subset subset expression if estimating PS only for the subset
##' @param ps_prefix string used as the prefix. Use to distinguish the re-estimated PS from the original PS.
##'
##' @return data_frame containing three additional columns for PS.
##'
##' @export
add_gps <- function(data,
                    formula,
                    family = multinomial(parallel = FALSE),
                    subset,
                    ps_prefix = "ps_") {
    assert_that("data.frame" %in% class(data))
    assert_that(class(formula) == "formula")
    assert_that(is.character(ps_prefix))
    assert_that(length(ps_prefix) == 1)

    ## Logical
    if (missing(subset)) {
        ## If NULL, using all the sample.
        subset_logical <- rep(TRUE, nrow(data))
    } else {
        ## Otherwise, subset appropriately.
        subset_logical <- eval(substitute(subset), envir = data)
    }

    ## Fit multinomial logistic regression
    res_vglm <- try(VGAM::vglm(formula = formula,
                               data    = data[subset_logical,],
                               family  = family))

    if (is.error(res_vglm)) {

        cat("## vglm errored. Checking for a design matrix rank deficiency.\n")
        ## Design matrix
        X <- model.matrix(object = formula, data = data)
        X_rank <- Matrix::rankMatrix(X)
        ## Look for columns that are linearly dependent.
        for (i in seq_len(ncol(X))) {
            if (X_rank == Matrix::rankMatrix(X[,-i])) {
                cat("##  Dropping ", colnames(X)[i], " does not reduce rank.\n")
            }
        }

        cat("## Returning data as is without adding PS..\n")
        return(data)

    } else {

        ## Calculate PS for all groups.
        ps_data <- as_data_frame(predict(res_vglm, newdata = data[subset_logical,], type = "response"))
        ## Add prefix
        names(ps_data) <- paste0(ps_prefix, names(ps_data))
        ##
        ps_data_full_size <- data_frame(v1 = rep(as.numeric(NA), nrow(data)),
                                        v2 = v1,
                                        v3 = v1)
        names(ps_data_full_size) <- names(ps_data)
        ## Put the calculated PS in appropriate rows
        ps_data_full_size[subset_logical,] <- ps_data
        ## Add to original data frame
        return(bind_cols(data, ps_data_full_size))
    }
}


###
### Functions for PS trimming
################################################################################


###
### Functions for PS weighting
################################################################################


###
### Functions for PS matching
################################################################################
