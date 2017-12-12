###
### Data analysis functions for three-group trimming paper
################################################################################

###
### Outcome analysis function
################################################################################

###   coxph analysis function (unadj, iptw, mw, ow) for anti-diabetics data
analyze_outcome_glm_log_linear <- function(data, formula, family) {

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
                                subset = (keep == 1),
                                robust = TRUE)),
                mw1 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = mw1,
                              subset = (keep == 1),
                              robust = TRUE)),
                mw2 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = mw2,
                              subset = (keep == 1),
                              robust = TRUE)),
                ow1 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = ow1,
                              subset = (keep == 1),
                              robust = TRUE)),
                ow2 = try(glm(formula = formula,
                              family = family,
                              data = data,
                              weights = ow2,
                              subset = (keep == 1),
                              robust = TRUE))
                )

    data_frame(adjustment = names(lst),
               model = lst)
}


###   Calculate estimate, SE, CI for coxph object (robust se is used if available via vcov())
calc_coxph_est_ci <- function(coxph_obj) {
    coefs <- coef(coxph_obj)
    ses <- sqrt(diag(vcov(coxph_obj)))
    data_frame(contrast = names(coefs),
               coef = coefs,
               se = ses,
               lower = coefs - qnorm(p = 0.975) * ses,
               upper = coefs + qnorm(p = 0.975) * ses)
}
