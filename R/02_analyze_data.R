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

###   glm outcome analysis function.
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


###   Calculate estimate, SE, CI for coxph object (robust se is used if available via vcov())
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
