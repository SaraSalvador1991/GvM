# Combine the rho-delta region plot with eight corresponding GvM plots.
# Run from the GvM folder: source("compare_gvm_modality.R")
source("plot_rho_delta.R")
source("plot_gvm.R")

# Compute the region grid once and reuse it for both output formats.
regions <- plot_rho_delta(draw = FALSE)
points_table <- data.frame(
  label = LETTERS[1:8],
  rho = c(rep(0.20, 3), rep(0.70, 3), rep(1.20, 2)),
  delta = c(0.40, 1.50, 2.60, 0.65, 1.57, 2.50, 0.85, 2.25)
)
points_table$bimodal <- mapply(is_bimodal, points_table$rho, points_table$delta)

plot_comparison <- function() {
  old <- par(no.readonly = TRUE)
  on.exit({layout(1); par(old)}, add = TRUE)
  layout(rbind(rep(1, 4), 2:5, 6:9), heights = c(1.8, 1, 1))
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
      mu1 = pi, mu2 = (pi - p$delta) %% pi,
      kappa1 = 4*p$rho*0.4, kappa2 = 0.4,
      main = paste(p$label, " | ", if (p$bimodal) "Bimodal" else "Unimodal"),
      interval = c(-pi/2, 3*pi/2)
    )
    mtext(sprintf("rho = %.2f; delta = %.2f", p$rho, p$delta),
          side = 3, line = 0.25, cex = 0.7)
  }
  invisible(densities)
}

png("compare_gvm_modality.png", width = 2400, height = 2200, res = 200)
densities <- plot_comparison()
dev.off()
pdf("compare_gvm_modality.pdf", width = 12, height = 11)
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
