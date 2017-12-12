###
### Data analysis functions for three-group trimming paper
################################################################################

###
### Outcome analysis function
################################################################################

###   coxph analysis function (unadj, iptw, mw, ow) for anti-diabetics data
##' Analyze outcome using glm
##'
##' .. content for details ..
##'
##' @param data data_frame assumed to be trimmed already.
##' @param formula
##' @param family
##'
##' @return data_frame with an adjustment method column and a model fit list column.
##'
##' @export
analyze_outcome_glm <- function(data, formula, family) {

    ## model = FALSE
    ## a logical value indicating whether _model frame_ should be kept.
    ## Drop the model frame to save space.

    lst <- list(unadj = try(glm(formula = formula,
                                family = family,
                                data = data,
                                model = FALSE)),
                iptw1 = try(glm(formula = formula,
                                family = family,
                                data = data,
                                weights = iptw1,
                                robust = TRUE,
                                model = FALSE)),
                iptw2 = try(glm(formula = formula,
                                family = family,
                                data = data,
                                weights = iptw2,
                                model = FALSE)),
                mw1 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = mw1,
                              model = FALSE)),
                mw2 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = mw2,
                              model = FALSE)),
                ow1 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = ow1,
                              model = FALSE)),
                ow2 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = ow2,
                              model = FALSE))
                )

    data_frame(adjustment = names(lst),
               ## List column
               model = lst)
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
