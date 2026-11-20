# GvM normalising constant G0, using the Fourier-Bessel series.
# delta = (mu1 - mu2) modulo pi; kappa1 and kappa2 are concentrations.
G0 <- function(delta, kappa1, kappa2, tolerance = 1e-12) {
  stopifnot(length(delta) == 1L, is.finite(delta),
            length(kappa1) == 1L, is.finite(kappa1), kappa1 >= 0,
            length(kappa2) == 1L, is.finite(kappa2), kappa2 >= 0,
            length(tolerance) == 1L, is.finite(tolerance), tolerance > 0)
  total <- besselI(kappa1, 0) * besselI(kappa2, 0)
  for (j in seq_len(10000L)) {
    coefficient <- 2 * besselI(kappa1, 2 * j) * besselI(kappa2, j)
    total <- total + coefficient * cos(2 * j * delta)
    if (!is.finite(total)) stop("Bessel series overflow: concentrations are too large.")
    # A negative or nearly zero cosine term does not imply convergence.
    # Stop using the unsigned coefficient instead of the signed summand.
    if (coefficient <= tolerance * abs(total)) return(total)
  }
  stop("G0 series did not converge within 10000 terms.")
}
