# Plot a GvM density for four parameters (angles in radians).
# Usage: source("plot_gvm.R"); plot_gvm(mu1, mu2, kappa1, kappa2)
source("gvm_normalizing_constant.R")

plot_gvm <- function(mu1, mu2, kappa1, kappa2, main = "",
                     interval = c(0, 2*pi)) {
  normalizer <- 2*pi*G0((mu1 - mu2) %% pi, kappa1, kappa2)
  density <- function(x) {
    exp(kappa1*cos(x - mu1) + kappa2*cos(2*(x - mu2))) / normalizer
  }
  theta <- seq(interval[1], interval[2], length.out = 2000)
  y <- density(theta)
  plot(theta, y, type = "l", lwd = 2, col = "#234E70", xaxt = "n",
       xlab = expression(theta), ylab = "Density", main = main,
       ylim = c(0, max(y)*1.08))
  if (identical(interval, c(0, 2*pi))) {
    axis(1, at = c(0, pi, 2*pi), labels = expression(0, pi, 2*pi))
  } else {
    axis(1, at = c(-pi/2, 0, pi, 3*pi/2),
         labels = expression(-pi/2, 0, pi, 3*pi/2))
  }
  invisible(density)
}

if (sys.nframe() == 0L) {
  plot_gvm(mu1 = pi/6, mu2 = pi/6 + pi/2, kappa1 = 1, kappa2 = 2)
}
