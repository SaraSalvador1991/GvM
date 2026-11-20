# Draw the GvM modality regions and return the computed grid.
# Usage: source("plot_rho_delta.R"); regions <- plot_rho_delta()

is_bimodal <- function(rho, delta, tolerance = 1e-7) {
  if (abs(cos(2 * delta)) < 1e-12 && rho < 0.5) return(TRUE)
  s <- sin(delta) * cos(delta)
  # polyroot() takes polynomial coefficients in ascending order.
  roots <- polyroot(c(s^2, 2*rho*s, rho^2 - 1, -4*rho*s, 1))
  all(abs(Im(roots)) < tolerance & abs(Re(roots)) <= 1 + tolerance)
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
