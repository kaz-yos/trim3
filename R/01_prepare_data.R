###
### Data preparation scripts for three-group CER
################################################################################

###
### Functions for PS estimation
################################################################################

###  Generalized PS Function
##' Add generalized PS
##'
##' .. content for details ..
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
                    ps_prefix) {
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

###   None (a dummy function to return 1's)
##' Dummy function to conduct no trimming at all.
##'
##' .. content for details ..
##'
##' @param A Treatment indicator vector having three levels indicated in \code{levels}.
##' @param ps0 PS vector for the first level.
##' @param ps1 PS vector for the second level.
##' @param ps2 PS vector for the third level.
##' @param levels Character vector holding three elements corresponding to levels in \code{A}
##' @param thres Threshold for the trimming strategy.
##'
##' @return A numeric indicator vector. 1 for kept and 0 for dropped.
##'
##' @export
trim_none <- function(A, ps0, ps1, ps2, levels, thres) {
    ## All of them has to be in [thres, 1.0].
    ## Treatment vector is not used.
    as.numeric(rep(1, length(ps0)))
}

###   Crump
##' Crump trimming function
##'
##' .. content for details ..
##'
##' @inheritParams trim_none
##'
##' @return A numeric indicator vector. 1 for kept and 0 for dropped.
##'
##' @export
trim_crump <- function(A, ps0, ps1, ps2, levels, thres) {
    ## All of them has to be in [thres, 1.0].
    ## Treatment vector is not used.
    as.numeric(ps0 >= thres & ps1 >= thres & ps2 >= thres)
}

###   Sturmer
##' Sturmer trimming function
##'
##' .. content for details ..
##'
##' @inheritParams trim_none
##'
##' @return A numeric indicator vector. 1 for kept and 0 for dropped.
##'
##' @export
trim_sturmer <- function(A, ps0, ps1, ps2, levels, thres) {

    ## Calculate the lower bounds as thres-quantile in respective group.
    l0 <- quantile(ps0[A == levels[1]], prob = thres)
    l1 <- quantile(ps1[A == levels[2]], prob = thres)
    l2 <- quantile(ps2[A == levels[3]], prob = thres)

    ## All three must be satisfied
    result <- as.numeric(ps0 >= l0 & ps1 >= l1 & ps2 >= l2)

    attributes(result) <- list(bounds = c(l0, l1, l2))

    result
}

###   Walker
##' Walker trimming function
##'
##' .. content for details ..
##'
##' @inheritParams trim_none
##'
##' @return A numeric indicator vector. 1 for kept and 0 for dropped.
##'
##' @export
trim_walker <- function(A, ps0, ps1, ps2, levels, thres) {

    ## Prevalence of treatment
    p0 <- mean(A == levels[1])
    p1 <- mean(A == levels[2])
    p2 <- mean(A == levels[3])

    ## Preference scores
    pi0 <- (ps0 / p0) / ((ps0 / p0) + (ps1 / p1) + (ps2 / p2))
    pi1 <- (ps1 / p1) / ((ps0 / p0) + (ps1 / p1) + (ps2 / p2))
    pi2 <- (ps2 / p2) / ((ps0 / p0) + (ps1 / p1) + (ps2 / p2))

    ## All three must be satisfied
    as.numeric(pi0 >= thres & pi1 >= thres & pi2 >= thres)
}

###   Generic trimming function (returns data with keep vector enclosed in a list)
##' Dummy function to conduct no trimming at all.
##'
##' .. content for details ..
##'
##' @param data data_frame
##' @param trim_method_name Function used for trimming that takes \code{A}, \code{ps0}, \code{ps1}, \code{ps2}, \code{levels}, and \code{thres} as arguments. Either one of the string: \code{\link{none}}, \code{\link{crump}}, \code{\link{sturmer}}, \code{\link{walker}}.
##' @param thres Threshold for the trimming strategy.
##' @param A_name Treatment indicator name in \code{data}.
##' @param ps0_name Name of the column in \code{data} for the PS for the first level.
##' @param ps1_name Name of the column in \code{data} for the PS for the second level.
##' @param ps2_name Name of the column in \code{data} for the PS for the third level.
##' @param levels Character vector holding three elements corresponding to levels in \code{A}
##'
##' @return trimmed data_frame with fewer rows
##'
##' @export
trim_data <- function(data, trim_method_name, thres,
                      A_name, ps0_name, ps1_name, ps2_name,
                      levels) {

    trim_fun <- get(paste0("trim_", trim_method))

    data$keep <- trim_fun(A   = unlist(data[,A_name]),
                          ps0 = unlist(data[,ps0_name]),
                          ps1 = unlist(data[,ps1_name]),
                          ps2 = unlist(data[,ps2_name]),
                          levels = levels,
                          thres = thres)

    data <- data %>%
        filter(keep == 1)

    ## It was a singleton list with data_frame for the empirical analyses.
    ## list(data)
    data
}


###
### Functions for PS weighting
################################################################################

###   Return one type of weight
##' Calculate one type of PS weight given PS
##'
##' .. content for details ..
##'
##' @inheritParams trim_none
##' @param weight_type Either one of \code{iptw}, \code{mw}, and \code{ow}.
##'
##' @return A numeric vector of individual PS weights.
##'
##' @export
calculate_weight <- function(A, ps0, ps1, ps2, levels, weight_type) {

    iptw <- (1/ps0 * (A == levels[1])) + (1/ps1 * (A == levels[2])) + (1/ps2 * (A == levels[3]))

    if (weight_type == "iptw") {
        ## IPTW
        return(iptw)

    } else if (weight_type == "mw") {
        ## Matching weights (MW)
        return(pmin(ps0, ps1, ps2) * iptw)

    } else if (weight_type == "ow") {
        ## Overlap weights (OW)
        return((ps0 * ps1 * ps2) * iptw)

    } else {
        stop("Invalid option for weight_type: ", weight_type)

    }
}

###   Add all weights
##' Calculate one type of PS weight given PS
##'
##' .. content for details ..
##'
##' @inheritParams trim_none
##' @param data data_frame
##' @param ps1_prefix Prefix for the PS estimated in the entire cohort.
##' @param ps2_prefix Prefix for the PS estimated in the trimmed cohort.
##'
##' @return data_frame containing IPTW, MW, and OW estimated in the entire cohort as well as the trimmed cohort.
##'
##' @export
add_all_weights <- function(data, A_name, levels, ps1_prefix = "ps1_", ps2_prefix = "ps2_") {

    ps1_names <- paste0("ps1_", levels)

    data$iptw1 <- NA
    data[data$keep == 1,]$iptw1 <- calculate_weight(A = unlist(data[data$keep == 1, A_name]),
                                                    ps0 = unlist(data[data$keep == 1, ps1_names[1]]),
                                                    ps1 = unlist(data[data$keep == 1, ps1_names[2]]),
                                                    ps2 = unlist(data[data$keep == 1, ps1_names[3]]),
                                                    levels = levels,
                                                    weight_type = "iptw")
    data$mw1 <- NA
    data[data$keep == 1,]$mw1 <- calculate_weight(A = unlist(data[data$keep == 1, A_name]),
                                                  ps0 = unlist(data[data$keep == 1, ps1_names[1]]),
                                                  ps1 = unlist(data[data$keep == 1, ps1_names[2]]),
                                                  ps2 = unlist(data[data$keep == 1, ps1_names[3]]),
                                                  levels = levels,
                                                  weight_type = "mw")
    data$ow1 <- NA
    data[data$keep == 1,]$ow1 <- calculate_weight(A = unlist(data[data$keep == 1, A_name]),
                                                  ps0 = unlist(data[data$keep == 1, ps1_names[1]]),
                                                  ps1 = unlist(data[data$keep == 1, ps1_names[2]]),
                                                  ps2 = unlist(data[data$keep == 1, ps1_names[3]]),
                                                  levels = levels,
                                                  weight_type = "ow")

    ps2_names <- paste0("ps2_", levels)
    ## Proceed if ps2_* are all available
    if (all(ps2_names %in% names(data))) {

        data$iptw2 <- NA
        data[data$keep == 1,]$iptw2 <- calculate_weight(A = unlist(data[data$keep == 1, A_name]),
                                                        ps0 = unlist(data[data$keep == 1, ps2_names[1]]),
                                                        ps1 = unlist(data[data$keep == 1, ps2_names[2]]),
                                                        ps2 = unlist(data[data$keep == 1, ps2_names[3]]),
                                                        levels = levels,
                                                        weight_type = "iptw")
        data$mw2 <- NA
        data[data$keep == 1,]$mw2 <- calculate_weight(A = unlist(data[data$keep == 1, A_name]),
                                                      ps0 = unlist(data[data$keep == 1, ps2_names[1]]),
                                                      ps1 = unlist(data[data$keep == 1, ps2_names[2]]),
                                                      ps2 = unlist(data[data$keep == 1, ps2_names[3]]),
                                                      levels = levels,
                                                      weight_type = "mw")
        data$ow2 <- NA
        data[data$keep == 1,]$ow2 <- calculate_weight(A = unlist(data[data$keep == 1, A_name]),
                                                      ps0 = unlist(data[data$keep == 1, ps2_names[1]]),
                                                      ps1 = unlist(data[data$keep == 1, ps2_names[2]]),
                                                      ps2 = unlist(data[data$keep == 1, ps2_names[3]]),
                                                      levels = levels,
                                                      weight_type = "ow")

    }
    ## Return data in the end
    return(data)
}


###
### Functions for PS matching
################################################################################


###
### Functions for coordinating the entire data preparation
################################################################################
##' One-stop function for data preparation
##'
##' Takes a data_frame and
##'
##' @param data data_frame
##' @param formula1 PS model formula for the entire cohort estimation.
##' @param family PS model family
##' @param ps_prefix1 PS variable prefix for the entire cohort estimation
##' @param ps_prefix2 PS variable prefix for the trimmed cohort estimation
##' @param df_trim_thres data_frame with columns \code{trim_method_name} and \code{thres} to indicate
##' @param A_name
##' @param ps0_name
##' @param ps1_name
##' @param ps2_name
##' @param levels
##' @param data_frame for one iteration
##'
##' @return nested data_frame containing
##'
##' @export
prepare_data <- function(data,
                         formula1,
                         family,
                         ps_prefix1,
                         ps_prefix2,
                         df_trim_thres = tribble(
                             ~trim_method_name, ~thres,
                             "none", 0,
                             "crump", 0.1,
                             "crump", 0.07,
                             "crump", 0.03,
                             "sturmer", 0.05,
                             "sturmer", 0.033,
                             "sturmer", 0.015,
                             "walker", 0.30,
                             "walker", 0.20,
                             "walker", 0.10),
                         A_name,
                         levels) {

    ## Add first-step GPS based on the entire cohort.
    data <- add_gps(data = data,
                    formula = formula1,
                    family = multinomial(parallel = FALSE),
                    ps_prefix = ps_prefix1)

    ## Trimming method and threshold configuration
    df_trim_thres <- tribble(
        ~trim_method_name, ~thres,
        "none", 0,
        "crump", 0.1,
        "crump", 0.07,
        "crump", 0.03,
        "sturmer", 0.05,
        "sturmer", 0.033,
        "sturmer", 0.015,
        "walker", 0.30,
        "walker", 0.20,
        "walker", 0.10)

    ## Trim data by various methods and hold in a nested data frame
    nested_df <- df_trim_thres %>%
        group_by(trim_method_name, thres) %>%
        mutate(trimmed_data = trim_data(data = data,
                                        ## String is used to retrieve the right function.
                                        trim_method_name = trim_method_name,
                                        thres = thres,
                                        A_name = A_name,
                                        ps0_name = paste0(ps_prefix1, levels)[1],
                                        ps1_name = paste0(ps_prefix1, levels)[2],
                                        ps2_name = paste0(ps_prefix1, levels)[3],
                                        levels = levels))

    ## Conduct PS re-estimation within each trimmed cohort.
    nested_df <- nested_df %>%
        group_by(trim_method_name, thres) %>%
        mutate(trimmed_data = map(trimmed_data, function(data) {
            add_gps(data = data,
                    formula = formula2,
                    family = multinomial(parallel = FALSE),
                    ps_prefix = ps_prefix2)
        }))

    ## Add weights
    nested_df <- nested_df %>%
        group_by(trim_method_name, thres) %>%
        mutate(trimmed_data = map(trimmed_data, function(data) {
            add_all_weights(data = data,
                            A_name = A_name,
                            levels = levels,
                            ps1_prefix = ps1_prefix,
                            ps2_prefix = ps2_prefix)
        }))

}
