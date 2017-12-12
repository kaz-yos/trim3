###
### Data analysis functions for three-group trimming paper
################################################################################

###
### Outcome analysis function
################################################################################

###   coxph analysis function (unadj, iptw, mw, ow) for anti-diabetics data
analyze_outcome_glm <- function(data, formula, family) {

    lst <- list(unadj = try(glm(formula = formula,
                                family = family,
                                data = data,
                                subset = (keep == 1))),
                iptw1 = try(glm(formula = formula,
                                family = family,
                                data = data,
                                weights = iptw1,
                                subset = (keep == 1),
                                robust = TRUE)),
                iptw2 = try(glm(formula = formula,
                                family = family,
                                data = data,
                                weights = iptw2,
                                subset = (keep == 1))),
                mw1 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = mw1,
                              subset = (keep == 1))),
                mw2 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = mw2,
                              subset = (keep == 1))),
                ow1 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = ow1,
                              subset = (keep == 1))),
                ow2 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = ow2,
                              subset = (keep == 1)))
                )

    data_frame(adjustment = names(lst),
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
