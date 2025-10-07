source("scripts/functions.R")

data_us <- read.csv("inputs/data/All/Mapping Police Violence_new.csv") %>%
  arrange(fips, year, month)

data_us[, 'tmean_lag1'] <- shift(data_us$tmean, 1, data_us$fips)
data_us[, 'prec_lag1'] <- shift(data_us$prec, 1, data_us$fips)

data_us[, 'tmean_lag0'] <- shift(data_us$tmean, 0, data_us$fips)
data_us[, 'prec_lag0'] <- shift(data_us$prec, 0, data_us$fips)

data_us[, 'tmean_lead1'] <- shift(data_us$tmean, -1, data_us$fips)
data_us[, 'prec_lead1'] <- shift(data_us$prec, -1, data_us$fips)
###########################################################
#The below code has been commented out because it
#takes ~15h to run. The stored data
#from thesese bootstrap runs will be loaded below.
#Uncomment out this section to run bootstraps and replace
#csv files with updated data.      
########################################################### 
# #define regression equation
fmla <- as.formula("rate_adj ~ tmean + tmean_lag1 + tmean_lead1 + prec + prec_lag1 + prec_lead1 | stateyear ")
feols_model <- feols(rate_adj ~ tmean + tmean_lag1 + tmean_lead1 +
                       prec + prec_lag1 + prec_lead1 | fipsmo + stateyear,
                     data = data_us,
                     cluster = ~fips)
n_clust <- length(unique(data_us$fips))
# bootstrap 1,000x with progress printing every 50 iterations
# est_us <- sapply(1:1000, function(x){
#   if(x %% 50 == 0) {
#     cat(sprintf("Progress: %d iterations completed (%.1f%%)\n", x, (x/1000)*100))
#   }
#   run_reg_bstrap(data_us, "fips")
# })
# write_csv(data.frame(est_us), path = "inputs/bootstrap_runs/BootstrapProjections_us_all.csv")

est_us <- read.csv("inputs/bootstrap_runs/BootstrapProjections_us_all.csv")

Tus <- read.csv("inputs/projections/TemperatureProjections_us.csv") %>%
  as_data_frame()

data_us <- read.csv("inputs/data/All/Mapping Police Violence_new.csv.csv")
data_us$yr <- as.numeric(as.character(data_us$year))
yrs_data <- 2013:2024
ind <- data_us$yr %in% yrs_data
data_us_sub <- data_us[ind,]
bus <- weighted.mean(data_us_sub$rate_adj, data_us_sub$popw, na.rm = TRUE)  

##############################################
### Calculate the cumulative excess mortality from 2020 to 2050####
### Assumption: Temperature increases linearly, baseline rate remains unchanged####
##############################################
pp <- read_csv("inputs/projections/PopulationProjections.csv")
ph <- read_csv("inputs/projections/PopulationProjections_Historical.csv")

p_hist <- ph[ph[,3] == "United States of America", paste0(2000:2015)]
popUS_2000_2015 <- as.numeric(gsub(" ", "", unlist(p_hist))) * 1000  
p_proj <- pp[pp[,3] == "United States of America", paste0(2016:2050)]
popUS_2016_2050 <- as.numeric(gsub(" ", "", unlist(p_proj))) * 1000  
popUS_full <- c(popUS_2000_2015, popUS_2016_2050)  

startYear <- 2020
endYear <- 2050
yr <- startYear:endYear

popUS_prediction <- popUS_full[(which(2000:2050 == startYear)):(which(2000:2050 == endYear))]
if(length(popUS_prediction) != (endYear - startYear + 1)){
  stop("Population data length does not match the projection period.")
}

#Calculate the additional mortality rate in the United States: Based on bootstrap results and temperature change forecasts (*12 converted to annual data, per 100,000 people per year)
dUS <- as.matrix(est_us) %*% matrix(Tus$deltaT, nrow = 1) * 12

ex <- array(dim = c(dim(dUS), length(yr)))
for (y in seq_along(yr)) {
  z <- (yr[y] - 2000) / (2050 - 2000)  
  ex[,,y] <- dUS * z * popUS_prediction[y] / 100000
}

exs <- apply(ex, c(1, 2), sum)

bxu <- c(exs)
qu <- quantile(bxu, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
bxu[bxu < qu[1]] <- qu[1]
bxu[bxu > qu[3]] <- qu[3]
print(qu)

med <- apply(ex, 3, median, na.rm = TRUE)

ylim_lower <- -100  
ylim_upper <- 1300    

alpha_plot <- 0.03
cll <- apply(sapply("orange", col2rgb) / 255, 2, function(x) {
  rgb(x[1], x[2], x[3], alpha = alpha_plot)
})


tiff(
  filename = "outputs/raw_figures/Figure7.tiff",
  width = 6,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw" 
)
par(mar = c(5, 5.5, 4, 2.2))
plot(
  yr,
  cumsum(ex[1, 1, ]),
  type = "n",
  ylim = c(ylim_lower, ylim_upper),
  xlim = c(startYear, endYear + 2),
  las = 1,
  ylab = "Excess deaths of police violence (cumulative)",
  xlab = "Year",
  cex.axis = 1.05,
  cex.lab = 1.05,
  axes = FALSE,
  yaxs = "i"
)
total_iterations <- dim(ex)[1] * dim(ex)[2]
iteration_count <- 0
for (j in 1:dim(ex)[1]) {
  for (i in 1:dim(ex)[2]) {
    lines(yr, cumsum(ex[j, i, ]), col = cll, lwd = 1)
    iteration_count <- iteration_count + 1
    if (iteration_count %% 50 == 0) {
      cat("Completed", iteration_count, "of", total_iterations, "iterations (", round(iteration_count / total_iterations * 100, 2), "%)\n")
    }
  }
}

lines(yr, cumsum(med), col = "black", lwd = 3)
abline(h = 0, lty = 2)
boxplot(
  bxu,
  horizontal = FALSE,
  range = 0,
  at = endYear + 2,
  add = TRUE,
  col = "orange",
  axes = FALSE,
  boxwex = 2,
  lty = 1
)

axis(1, at = seq(startYear, endYear, 5), labels = FALSE, cex.axis = 1.2)
mtext(text = expression(italic("2020")), side = 1, at = 2020, line = 1, cex = 1.2, adj = 0.5, srt = 45)
mtext(text = expression(italic("2025")), side = 1, at = 2025, line = 1, cex = 1.2, adj = 0.5, srt = 45)
mtext(text = expression(italic("2030")), side = 1, at = 2030, line = 1, cex = 1.2, adj = 0.5, srt = 45)
mtext(text = expression(italic("2035")), side = 1, at = 2035, line = 1, cex = 1.2, adj = 0.5, srt = 45)
mtext(text = expression(italic("2040")), side = 1, at = 2040, line = 1, cex = 1.2, adj = 0.5, srt = 45)
mtext(text = expression(italic("2045")), side = 1, at = 2045, line = 1, cex = 1.2, adj = 0.5, srt = 45)
mtext(text = expression(italic("2050")), side = 1, at = 2050, line = 1, cex = 1.2, adj = 0.5, srt = 45)

axis(2,
     at = seq(ylim_lower, ylim_upper, 100),
     labels = seq(ylim_lower, ylim_upper, 100),
     las = 1,
     cex.axis = 1.2
)
dev.off()