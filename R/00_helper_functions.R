################################################################################
### Helper functions
##
## Created on: 2017-12-15
## Author: Kazuki Yoshida
################################################################################

###
### Simple math functions
################################################################################
expit <- function(x) {
    e_x <- exp(x)
    e_x / (1 + e_x)
}
##
logit <- function(p) {
    log(p / (1 - p))
}


###
### Error handling
################################################################################

###  Detect errors given by try()
## http://adv-r.had.co.nz/Exceptions-Debugging.html#condition-handling
##' Tests for objects of type \code{try-error}
##'
##' Tests for objects of type \code{try-error}
##'
##' @param x return object from a function wrapped by \code{try}.
##'
##' @return logical indicating whether the return object is an error object given by \code{try} or normal return object.
##'
##' @author Kazuki Yoshida
##'
##' @export
is.error <- function(x) {
    inherits(x, "try-error")
}
