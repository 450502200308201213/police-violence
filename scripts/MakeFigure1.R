source("scripts/functions.R")
data_us <- read.csv("inputs/data/Mapping Police Violence_new")

###########################################################
#The below code has been commented out because it
#takes ~30 minutes to run. The stored data
#from thesese bootstrap runs will be read-in below.
#Uncomment out this section to run bootstraps and replace
#file `BootstrapMainModel_us.csv' with updated data.      
###########################################################   
# run bootstrap of base model for US (county-month + state-year FE). this takes a while. 
# cc <- unique(data_us$fips)
# set.seed(1)
# xx <- -17:33
# spline_out <- NULL
# for (i in 1:1000) {
#   tryCatch({
#     samp <- data.frame(fips=sample(cc, length(cc), replace=TRUE))
#     subdata <- inner_join(data_us, samp, by="fips", relationship="many-to-many")
#     reg <- felm(rate_adj ~ ns(tmean, knots=c(0, 10, 20)) + prec | stateyear | 0 | fips,
#                 weights=subdata$popw, data=subdata)
#     spline_est <- as.numeric(t(as.matrix(coef(reg)[1:4])) %*%
#                                t(matrix(nrow=length(xx), ncol=4, data=ns(xx, knots=c(0, 10, 20)))))
#     spline_out <- rbind(spline_out, spline_est)
#   }, error=function(e){
#     print(paste("Error in iteration", i, ":", e))
#   })
#   print(i)
# }
# 
# if (!is.null(spline_out) && nrow(spline_out) > 0) {
#   out_df <- data.frame(t(spline_out))
#   colnames(out_df) <- paste0("iter_", 1:nrow(spline_out))
#   out_df$temp <- xx
#   write_csv(out_df, path="inputs/bootstrap_runs/BootstrapSplineModel.csv")
# } else {
#   stop("Bootstrap failed to generate any valid results.")
# }
# out_df <- read_csv("inputs/bootstrap_runs/BootstrapSplineModel.csv")

### Run alternative specifications ###
data_us$yearmonth <- data_us$yr * 100 + data_us$month
xx <- -17:33

mod1_feols <- feols(rate_adj ~ tmean + prec | stateyear+fipsmo,
                    data = data_us,
                    weights = data_us$popw,
                    cluster = ~ fips)
yy <- data.frame(xx, stateyear = coef(mod1_feols)[1] * xx)

# mod2
mod2_feols <- feols(rate_adj ~ poly(tmean, 3, raw = TRUE) + prec | stateyear+fipsmo,
                    data = data_us,
                    weights = data_us$popw,
                    cluster = ~ fips)
yy <- data.frame(yy, poly = as.numeric(t(as.matrix(coef(mod2_feols)[1:3])) %*% t(matrix(nrow = length(xx), ncol = 3, data = poly(xx, 3, raw = TRUE)))))

# mod3
mod3_feols <- feols(rate_adj ~ ns(tmean, knots = c(0, 10, 20)) + prec | stateyear+fipsmo,
                    data = data_us,
                    weights = data_us$popw,
                    cluster = ~ fips)
yy <- data.frame(yy, spline = as.numeric(t(as.matrix(coef(mod3_feols)[1:4])) %*% t(matrix(nrow = length(xx), ncol = 4, data = ns(xx, knots = c(0, 10, 20))))))

# mod4
mod4_feols <- feols(rate_adj ~ bs(tmean, df = 4) + prec | stateyear+fipsmo,
                    data = data_us,
                    weights = data_us$popw,
                    cluster = ~ fips)
yy <- data.frame(yy, bspline = as.numeric(t(as.matrix(coef(mod4_feols)[1:4])) %*% t(matrix(nrow = length(xx), ncol = 4, data = bs(xx, df = 4)))))

# mod5
mod5_feols <- feols(rate_adj ~ ns(tmean, df = 8) + prec | stateyear+fipsmo,
                    data = data_us,
                    weights = data_us$popw,
                    cluster = ~ fips)
yy <- data.frame(yy, spline7 = as.numeric(t(as.matrix(coef(mod5_feols)[1:8])) %*% t(matrix(nrow = length(xx), ncol = 8, data = ns(xx, df = 8)))))

### Plot Figure ###
png(file="outputs/raw_figures/Figure1.png",
    height=8, width=6, units="in", res=600)
par(mfrow=c(1, 1),mar=c(5, 5, 4, 2) + 0.1)

#read-in bootstrapped runs from above to get CI
boot_spline <- read.csv("inputs/bootstrap_runs/BootstrapSplineModel.csv")
spline_est <- as.matrix(boot_spline[, -ncol(boot_spline)])
spline_est <- t(spline_est)
spline_est <- spline_est - spline_est[, which(xx==10)]
ci_spline <- apply(spline_est, 2, function(x) quantile(x, probs=c(0.025, 0.975))) / weighted.mean(data_us$rate_adj, data_us$popw) * 100

yyr <- yy[, c("xx", "poly", "spline", "spline7", "bspline")]
nn <- dim(yyr)[2]
br <- weighted.mean(data_us$rate_adj, data_us$popw)

# Check for invalid br
if (is.na(br) || br == 0) {
  stop("Weighted mean (br) is NA or 0, which will cause division issues.")
}

for (i in 2:nn) {
  yyr[, i] <- (yyr[, i] - yyr[yyr$xx==10, i]) / br * 100
}

# Check for NA or Inf in yyr
if (any(is.na(yyr)) || any(is.infinite(as.matrix(yyr)))) {
  print("NA or Inf values found in yyr:")
  print(which(is.na(yyr) | is.infinite(as.matrix(yyr)), arr.ind=TRUE))
  stop("Please fix NA or Inf values in yyr before proceeding.")
}

adjusted_ci_spline <- ci_spline
for (i in 1:length(xx)) {
  spline_value <- yyr$spline[i]
  lower_ci <- ci_spline[1, i]
  upper_ci <- ci_spline[2, i]
  ci_range <- (upper_ci - lower_ci) / 2.0
  adjusted_ci_spline[1, i] <- spline_value - ci_range
  adjusted_ci_spline[2, i] <- spline_value + ci_range
}

ylim_lower <- min(c(adjusted_ci_spline[1, ], as.numeric(as.matrix(yyr[, 2:nn]))), na.rm=TRUE)
ylim_upper <- max(c(adjusted_ci_spline[2, ], as.numeric(as.matrix(yyr[, 2:nn]))), na.rm=TRUE)

plot(1, type="n", las=1, ylim=c(ylim_lower, ylim_upper), xlim=c(-20, 40),
     xlab="Monthly average temperature(°C)",
     ylab="Percentage change in death rate of police violence (%)",
     cex.lab=1.2
)
polygon(c(xx, rev(xx)), c(adjusted_ci_spline[1, ], rev(adjusted_ci_spline[2, ])), col="lightblue", border=NA)
abline(h=0, col="red", lty=3, lwd=1)

clz <- c("red", "blue", "purple", "green") 
lp <- c(1, 1, 1, 1) 

for (i in 2:nn) {
  print(paste("Drawing curve", i-1, ":", names(yyr)[i], "with color", clz[i-1], "and lty", lp[i-1]))
  lines(xx, yyr[, i], col=clz[i-1], lty=lp[i-1], lwd=2)
}

bb <- 50
ht <- hist(data_us$tmean, breaks=bb, plot=F)
bb <- length(ht$breaks)
rect(ht$breaks[1:(bb-1)], ylim_lower, ht$breaks[2:bb], ylim_lower + ht$counts/max(ht$counts)*35, col="deepskyblue3")

legend("topleft", legend=c("Confidence Interval(spline)", "poly", "spline", "spline7", "bspline"),
       fill=c("lightblue", "red", "blue", "purple", "green"),
       border=NA, bty="n", cex=1.38)

dev.off()