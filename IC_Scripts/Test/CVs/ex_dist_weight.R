# Distances from 0 to 100
d <- seq(0, 100, by = 1)

# Try different sigma values
sigmas <- c(10, 20, 30, 40)

# Plot the y=0 (no spill) line (always 1)
plot(d, rep(1, length(d)), type = "l", lwd = 2, col = "red", lty = 2,
     ylim = c(0, 1), xlab = "Distance d", ylab = "Weight v(d)",
     main = "Distance-based weights for spill vs no spill")

# Add curves for y=1 (spill) case, one for each sigma
cols <- c("blue", "darkgreen", "purple", "orange")
for (j in seq_along(sigmas)) {
  sigma <- sigmas[j]
  v_spill <- exp(-d / sigma)   # spill case (y = 1)
  lines(d, v_spill, lwd = 2, col = cols[j])
}

legend("topright",
       legend = c("y = 0", paste("y = 1, sigma =", sigmas)),
       col = c("red", cols),
       lwd = 2, lty = c(2, rep(1, length(sigmas))))
