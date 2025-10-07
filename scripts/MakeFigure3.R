source("scripts/functions.R")

# Load US data and interact with time effects
data_us <- read.csv("inputs/data/all/Mapping Police Violence_new.csv")
#data_us <- read.csv("inputs/data/White/White_new.csv")
#data_us <- read.csv("inputs/data/Unknown/Unknown race_new.csv")
#data_us <- read.csv("inputs/data/Native Hawaiian and Pacific Islander/Native Hawaiian and Pacific Isl_new.csv")
#data_us <- read.csv("inputs/data/Hispanic/Hispanic_new.csv")
#data_us <- read.csv("inputs/data/Black/Black_new.csv")
#data_us <- read.csv("inputs/data/Asian/Asian_new.csv")
#data_us <- read.csv("inputs/data/American Indian and Alaska Native/American Indian and Alaska Nati_new.csv")

yrs <- 2013:2024
data_us$year <- as.factor(data_us$year)

# Run model using feols
mod <- feols(rate_adj ~ tmean*year + prec*year | stateyear+fipsmo, data = data_us, weights = ~popw, cluster = ~fips)
coef <- summary(mod)$coefficients
vars <- grep("tmean", names(coef))

# Pull out and process results  
vcov <- vcov(mod)[vars, vars]
coef <- coef[vars]
b <- coef[1]
se <- summary(mod)$se[vars][1]

for (k in 2:length(vars)) {
  b <- c(b, coef[k] + coef[1])
  covv <- vcov[c(1, k), c(1, k)]
  se <- c(se, sqrt(sum(covv)))
}

cilo <- b - 1.96 * se
cihi <- b + 1.96 * se

res <- data.frame(yr = yrs, b, se, cilo, cihi)
avg <- data_us %>% group_by(yr) %>% summarize(avg = weighted.mean(rate_adj, popw))
res <- left_join(res, avg, by = "yr")

# Run simpler model for effect size
mod <- feols(rate_adj ~ tmean + prec | stateyear+fipsmo, data = data_us, weights = ~popw, cluster = ~fips)
br <- weighted.mean(data_us$rate_adj, data_us$popw)
eff <- coef(mod)["tmean"] / br * 100

# Plot panels
pdf(file = "outputs/raw_figures/Figure3.pdf", height = 5, width = 6, useDingbats = FALSE)
par(mfrow = c(1, 1))

plot(1, type = "n", las = 1, ylim = c(-20, 20), xlim = c(2013, 2024), ylab = "Effect on death rate of police violence (%)", xlab = "Year",yaxt = "n")
axis(2,at = seq(-20, 20, 5),labels = seq(-20, 20, 5),las = 1,cex.axis = 0.8)
abline(h = 0, col = "red", lty = 2)
segments(res$yr, res$cilo / res$avg * 100, res$yr, res$cihi / res$avg * 100, lwd = 2, col = "grey")
points(res$yr, res$b / res$avg * 100, pch = 19, cex = 1)

dev.off()