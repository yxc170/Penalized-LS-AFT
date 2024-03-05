#' Function to calculate the 1st and the 2nd moments of the residuals
#'
#' This function is copied from the aftgee package. Not intend to export.
#'
#' @importFrom survival survfit
#'
#' @keywords internal
#' @noRd
eRes <- function (e, delta, z = rep(1, length(e))) {
    nobs <- length(e)
    ord <- order(e)
    ei <- e[ord]
    deltai <- delta[ord]
    zi <- z[ord]
    dummy <- 1:nobs
    tmp <- survfit(Surv(ei, deltai) ~ 1, weights = zi)
    Shat <- with(tmp, approx(time, surv, ei))$y
    edif <- c(diff(ei), 0)
    ehat <- rev(cumsum(rev(edif * Shat)))
    inpt <- mean(ehat)
    ehat2 <- rev(cumsum(rev(ei * edif * Shat)))
    ehat <- ehat/Shat + ei
    ehat2 <- 2 * ehat2/Shat + ei^2
    ehat[is.na(ehat)] <- ei[is.na(ehat)]
    ehat2[is.na(ehat2)] <- ei[is.na(ehat2)]^2
    ehat2[which(ehat2 < 0)] <- NaN
    eres <- ehat
    eres2 <- ehat2
    eres[dummy[ord]] <- ehat
    eres2[dummy[ord]] <- ehat2
    return(list(eres, eres2, inpt))
}
