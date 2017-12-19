###
### Data analysis functions for three-group trimming paper
################################################################################

###
### Outcome analysis function
################################################################################

coef_and_vcov <- function(model, vcov_fun) {

    coef_vector <- coef(model)
    vcov_matrix <- vcov_fun(model)

    list(coef = coef_vector,
         vcov = vcov_matrix)
}

###  Data augmentation for counterfactual analyses
##' Create an augmented dataset of counterfactuals
##'
##' .. content for details ..
##'
##' @param data data_frame containing outcome, treatment, counterfactual outcomes under each one of three treatment values (0,1,2).
##' @param outcome_name name of the outcome variable.
##' @param counter_names names of the three counterfactual outcome variables
##' @param A_name name of the treatment variable
##' @param A_levels three levels for the treatment variable
##' @param ps1_prefix Prefix for the PS estimated in the entire cohort
##' @param ps2_prefix Prefix for the PS estimated in the trimmed cohort
##'
##' @return data_frame that is three times larger containing counterfactual outcomes for each individuals, weights are calculated from the true PS.
##'
##' @export
augment_counterfactuals <- function(data, outcome_name, counter_names, A_name, A_levels, ps1_prefix, ps2_prefix) {

    ## Clone
    data0 <- data
    data1 <- data
    data2 <- data
    ## Assign counterfactual outcomes
    data0[,outcome_name] <- data0[,counter_names[1]]
    data1[,outcome_name] <- data1[,counter_names[2]]
    data2[,outcome_name] <- data2[,counter_names[3]]
    ## Assign treatment
    data0[,A_name] <- A_levels[1]
    data1[,A_name] <- A_levels[2]
    data2[,A_name] <- A_levels[3]
    ## Combine into one
    data_aug <- bind_rows(data0,
                          data1,
                          data2)
    ## Replace first-stage PS with true PS
    ps1_names <- paste0(ps1_prefix, A_levels)
    data[, ps1_names[1]] <- data[, "pA0"]
    data[, ps1_names[2]] <- data[, "pA1"]
    data[, ps1_names[3]] <- data[, "pA2"]
    ## No need for re-estimated PS as they are true PS.
    ## Calculate weights from the true PS
    data$iptw1 <- calculate_weight(A = unlist(data[, A_name]),
                                   ps0 = unlist(data[, ps1_names[1]]),
                                   ps1 = unlist(data[, ps1_names[2]]),
                                   ps2 = unlist(data[, ps1_names[3]]),
                                   levels = A_levels,
                                   weight_type = "iptw")
    data$iptw2 <- data$iptw1
    data$mw1 <- calculate_weight(A = unlist(data[, A_name]),
                                 ps0 = unlist(data[, ps1_names[1]]),
                                 ps1 = unlist(data[, ps1_names[2]]),
                                 ps2 = unlist(data[, ps1_names[3]]),
                                 levels = A_levels,
                                 weight_type = "mw")
    data$mw2 <- data$mw1
    data$ow1 <- calculate_weight(A = unlist(data[, A_name]),
                                 ps0 = unlist(data[, ps1_names[1]]),
                                 ps1 = unlist(data[, ps1_names[2]]),
                                 ps2 = unlist(data[, ps1_names[3]]),
                                 levels = A_levels,
                                 weight_type = "ow")
    data$ow2 <- data$ow1
    ## Return the completed data
    data_aug
}


###  glm outcome analysis function.
##' Analyze outcome using glm
##'
##' .. content for details ..
##'
##' @param data data_frame which is assumed to have been trimmed already. It has to contain weight variables (\code{iptw1}, \code{iptw2}, \code{mw1}, \code{mw2}, \code{ow1}, and \code{ow2})
##' @param formula outcome analysis formula. It should only contain the outcome and the exposure of interest.
##' @param family \code{glm} family statement.
##' @param data_aug augmented data including counterfactual clones of each individual.
##'
##' @return data_frame with an adjustment method column and list columns for the coefficient and variance-covariance matrix (robust one for weighted data).
##'
##' @export
analyze_outcome_glm <- function(data, formula, family, data_aug) {

    lst <- list(unadj = try(coef_and_vcov(glm(formula = formula,
                                              family = family,
                                              data = data),
                                          vcov_fun = vcov)),
                iptw1 = try(coef_and_vcov(glm(formula = formula,
                                              family = family,
                                              data = data,
                                              weights = iptw1),
                                          vcov_fun = sandwich::sandwich)),
                iptw2 = try(coef_and_vcov(glm(formula = formula,
                                              family = family,
                                              data = data,
                                              weights = iptw2),
                                          vcov_fun = sandwich::sandwich)),
                mw1 = try(coef_and_vcov(glm(formula = formula,
                                            family = family,
                                            data = data,
                                            weights = mw1),
                                        vcov_fun = sandwich::sandwich)),
                mw2 = try(coef_and_vcov(glm(formula = formula,
                                            family = family,
                                            data = data,
                                            weights = mw2),
                                        vcov_fun = sandwich::sandwich)),
                ow1 = try(coef_and_vcov(glm(formula = formula,
                                            family = family,
                                            data = data,
                                            weights = ow1),
                                        vcov_fun = sandwich::sandwich)),
                ow2 = try(coef_and_vcov(glm(formula = formula,
                                            family = family,
                                            data = data,
                                            weights = ow2),
                                        vcov_fun = sandwich::sandwich))
                )

    ## Add counterfactual results if augmented datasets are available.
    if (!missing(data_aug)) {

        lst_aug <- list(unadj = try(coef_and_vcov(glm(formula = formula,
                                                      family = family,
                                                      data = data_aug),
                                                  vcov_fun = vcov)),
                        iptw1 = try(coef_and_vcov(glm(formula = formula,
                                                      family = family,
                                                      data = data_aug,
                                                      weights = iptw1),
                                                  vcov_fun = sandwich::sandwich)),
                        iptw2 = try(coef_and_vcov(glm(formula = formula,
                                                      family = family,
                                                      data = data_aug,
                                                      weights = iptw2),
                                                  vcov_fun = sandwich::sandwich)),
                        mw1 = try(coef_and_vcov(glm(formula = formula,
                                                    family = family,
                                                    data = data_aug,
                                                    weights = mw1),
                                                vcov_fun = sandwich::sandwich)),
                        mw2 = try(coef_and_vcov(glm(formula = formula,
                                                    family = family,
                                                    data = data_aug,
                                                    weights = mw2),
                                                vcov_fun = sandwich::sandwich)),
                        ow1 = try(coef_and_vcov(glm(formula = formula,
                                                    family = family,
                                                    data = data_aug,
                                                    weights = ow1),
                                                vcov_fun = sandwich::sandwich)),
                        ow2 = try(coef_and_vcov(glm(formula = formula,
                                                    family = family,
                                                    data = data_aug,
                                                    weights = ow2),
                                                vcov_fun = sandwich::sandwich))
                        )
        names(lst_aug) <- paste0(names(lst_aug), "_t")

        lst <- c(lst, lst_aug)
    }

    data_frame(adjustment = names(lst),
               ## List column
               coef = lapply(lst, `[[`, "coef"),
               vcov = lapply(lst, `[[`, "vcov"))
}


###  Calculate estimate, SE, CI for coxph object (robust se is used if available via vcov())
calc_est_ci <- function(model, vcov_fun, alpha) {
    ## Coefficients
    coefs <- coef(model)
    ## Standard errors using vcov_fun.
    ses <- sqrt(diag(vcov_fun(model)))
    data_frame(contrast = names(coefs),
               coef = coefs,
               se = ses,
               lower = coefs - qnorm(p = 1 - alpha / 2) * ses,
               upper = coefs + qnorm(p = 1 - alpha / 2) * ses)
}
