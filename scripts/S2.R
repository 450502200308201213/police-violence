race_lab  <- c("All",
               "White",
               "Unknown",
               "Native Hawaiian\nand Pacific Islander",
               "Hispanic",
               "Black",
               "Asian",
               "American Indian\nand Alaska Native")
est       <- c(1.873389, 1.229075, 5.658206, -3.006437,
               1.369409, 2.514795, -1.196291, 1.072319)
ci_lower  <- c(0.858543, -0.3353263, 2.079918, -16.02043,
               0.3046687, 0.437864, -8.550158, -8.3949)
ci_upper  <- c(2.888235, 2.793476, 9.236494, 10.00756,
               2.434149, 4.591727, 6.157576, 10.53954)
x_pos     <- seq(2013, 2024, length.out = 8)
#The data above are calculated based on the average effect values for different racial groups from MakeFigure3

######Plot######
pdf(file = "outputs/raw_figures/S2.pdf",
    height = 5, width = 10, useDingbats = FALSE)

par(mfrow = c(1, 1),
    mar = c(8, 5, 2, 2) + 0.1,   
    las = 1)                     

plot(1, type = "n",
     ylim = c(-15, 15),
     xlim = c(2013, 2024),
     ylab = "Effect on death rate of police violence (%)",
     xlab = "",
     xaxt = "n")

axis(1, at = x_pos, labels = race_lab, cex.axis = 0.75, las = 1)

abline(h = 0, col = "red", lty = 2)

segments(x_pos, pmax(ci_lower, -15),
         x_pos, pmin(ci_upper,  15),
         lwd = 2, col = "grey")
points(x_pos, est, pch = 19, cex = 1.1)
dev.off()