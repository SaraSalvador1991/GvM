# Draw the GvM modality regions and return the computed grid.
# Usage: source("plot_rho_delta.R"); regions <- plot_rho_delta()

# Analytic classification from Sections 2.1-2.2 of Salvador and Gatto (2021).
# Count proper maxima, not polynomial roots with their algebraic multiplicity.
# Points within numerical tolerance of a transition are treated as unimodal.
is_bimodal <- function(rho, delta, tolerance = 1e-10) {
  stopifnot(length(rho) == 1L, is.finite(rho), rho > 0,
            length(delta) == 1L, is.finite(delta),
            length(tolerance) == 1L, is.finite(tolerance), tolerance > 0)
  delta <- delta %% pi

  # delta = 0 or pi/2: the double root x = 0 must be mapped back to angles.
  if (abs(sin(2*delta)) <= tolerance) return(rho < 1 - tolerance)

  # delta = pi/4 or 3*pi/4: two double roots; one leaves [-1, 1]
  # above rho = 1/2. At equality, merged angles do not give a second maximum.
  if (abs(cos(2*delta)) <= tolerance) return(rho < 0.5 - tolerance)

  if (rho < 0.5) return(TRUE)
  if (rho >= 1) return(FALSE)

  # Sign of Delta_ft, without division by 27*rho^4:
  # 4*(1-rho^2)^3 - 27*rho^4*sin(2*delta)^2.
  # A zero corresponds to the generic transition: a stationary inflection,
  # rather than a second proper maximum. No quartic root solver is needed.
  left <- 4*(1 - rho^2)^3
  right <- 27*rho^4*sin(2*delta)^2
  left - right > tolerance * max(1, left, right)
}

plot_rho_delta <- function(regions = NULL, draw = TRUE, main = "") {
  if (is.null(regions)) {
    rho <- seq(0.005, 2, length.out = 300)
    delta <- seq(0.001, pi - 0.001, length.out = 400)
    regions <- list(rho = rho, delta = delta,
                    bimodal = outer(rho, delta, Vectorize(is_bimodal)))
  }
  if (draw) {
    image(regions$rho, regions$delta, regions$bimodal * 1,
          col = c("white", "grey75"), breaks = c(-0.5, 0.5, 1.5),
          xlab = expression(rho), ylab = expression(delta),
          xaxt = "n", yaxt = "n", main = main)
    axis(1, at = c(0.5, 1, 1.5, 2))
    axis(2, at = c(0, pi/4, pi/2, 3*pi/4, pi),
         labels = expression(0, pi/4, pi/2, 3*pi/4, pi))
  }
  invisible(regions)
}

# Running this file with Rscript draws the standalone example.
if (sys.nframe() == 0L) plot_rho_delta()
