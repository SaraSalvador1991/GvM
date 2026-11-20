# Combine the rho-delta region plot with nine corresponding GvM plots.
# Run from the GvM folder: source("compare_gvm_modality.R")
source("plot_rho_delta.R")
source("plot_gvm.R")

# Compute the region grid once and reuse it for both output formats.
# mu1 and kappa2 are explicit: rho and delta alone do not fix the density.
mu1 <- pi
kappa2 <- 1/6  # Estimated from the density scales in the reference figure.
regions <- plot_rho_delta(draw = FALSE)
points_table <- data.frame(
  label = LETTERS[1:9],
  rho = rep(c(0.3, 0.75, 1.5), each = 3),
  delta = rep(c(3*pi/4, pi/2, pi/4), 3),
  delta_label = rep(c("3*pi/4", "pi/2", "pi/4"), 3)
)

# Use higher concentrations for A and E to make their two modes easier to see.
# rho and delta remain identical to the reference figure.
points_table$kappa2 <- kappa2
points_table$kappa2[points_table$label == "A"] <- 0.75
points_table$kappa2[points_table$label == "E"] <- 1.5
points_table$bimodal <- mapply(is_bimodal, points_table$rho, points_table$delta)

plot_comparison <- function() {
  old <- par(no.readonly = TRUE)
  on.exit({layout(1); par(old)}, add = TRUE)
  layout(rbind(rep(1, 3), c(2, 5, 8), c(3, 6, 9), c(4, 7, 10)),
         heights = c(1.8, 1, 1, 1))
  par(mar = c(3.6, 4, 3.4, 1), mgp = c(2.3, 0.7, 0))
  plot_rho_delta(regions, main = "GvM: parameter regions and corresponding densities")
  points(points_table$rho, points_table$delta, pch = 20)
  text(points_table$rho, points_table$delta, points_table$label, pos = 4, font = 2)
  legend("topright", c("Bimodal", "Unimodal"), fill = c("grey75", "white"), bty = "n")

  densities <- vector("list", nrow(points_table))
  for (i in seq_len(nrow(points_table))) {
    par(mar = c(3.4, 3.6, 3.1, 0.8), cex.main = 0.95)
    p <- points_table[i, ]
    densities[[i]] <- plot_gvm(
      mu1 = mu1, mu2 = (mu1 - p$delta) %% pi,
      kappa1 = 4*p$rho*p$kappa2, kappa2 = p$kappa2,
      main = sprintf("Point %s = (%.2g, %s)", p$label, p$rho, p$delta_label),
      interval = c(0, 2*pi)
    )

  }
  invisible(densities)
}

png("compare_gvm_modality.png", width = 2200, height = 2700, res = 200)
densities <- plot_comparison()
dev.off()
pdf("compare_gvm_modality.pdf", width = 11, height = 13.5)
plot_comparison()
dev.off()
if (interactive()) plot_comparison()

# Check normalisation and the number of modes on a circular grid.
theta <- seq(0, 2*pi, length.out = 20001)[-20001]
count_modes <- function(f) {
  y <- f(theta)
  n <- length(y)
  sum(y > y[c(n, seq_len(n - 1))] & y > y[c(2:n, 1)])
}
stopifnot(all(vapply(densities, count_modes, integer(1)) ==
                ifelse(points_table$bimodal, 2L, 1L)))
for (f in densities) stopifnot(abs(integrate(f, 0, 2*pi)$value - 1) < 1e-7)

message("Comparison complete. PNG and PDF saved in: ", getwd())
