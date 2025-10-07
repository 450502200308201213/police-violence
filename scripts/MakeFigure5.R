source("scripts/functions.R")
# Load US data and filter for 2013–2015
data_us <- read.csv("inputs/data/All/Mapping Police Violence_new") %>%
  filter(year %in% 2013:2024) #The value can be adjusted within the specified range to match the desired year

###########################################################
# The below code has been commented out because it
# takes several hours to run. The stored data
# from these bootstrap runs will be loaded below.
# Uncomment this section to run bootstraps and replace
# file `BootstrapStateEffects_main_us.csv` and `BootstrapStateEffects_state_us.csv` with updated data.
###########################################################
# B <- 1000
# n1 <- nrow(data_us)
# 
# data_us$bstrap_id <- 1:nrow(data_us)
# 
# bstrap1_main <- rep(NA, B); names(bstrap1_main) <- paste("b", 1:B, sep = "")
# bstrap1_state <- matrix(nrow = length(unique(data_us$state)), ncol = B)
# rownames(bstrap1_state) <- sort(unique(data_us$state))
# colnames(bstrap1_state) <- paste("b", 1:B, sep = "")
# 
# us_states <- data.frame(state = sort(unique(data_us$state)))
# 
# # Bootstrap
# for (b in 1:B) {
#   data_us.sample <- left_join(data.frame(bstrap_id = sample(x = data_us$bstrap_id, size = n1, replace = TRUE)), data_us, "bstrap_id")
#   while(length(unique(data_us.sample$state)) < length(unique(data_us$state))) {
#     data_us.sample <- left_join(data.frame(bstrap_id = sample(x = data_us$bstrap_id, size = n1, replace = TRUE)), data_us, "bstrap_id")
#   }
#   
#   # US model using feols
#   mod1_main <- feols(rate_adj ~ tmean + prec | stateyear, data = data_us.sample, weights = ~popw)
#   mod1_state <- feols(rate_adj ~ tmean * as.factor(state) + prec * as.factor(state) | year, data = data_us.sample, weights = ~popw)
#   
#   bstrap1_main[b] <- coef(mod1_main)["tmean"]
#   coef <- coef(mod1_state)[grep("tmean", names(coef(mod1_state)))]
#   coef[2:length(coef)] <- coef[1] + coef[2:length(coef)]
#   bstrap1_state[, b] <- as.numeric(coef)
#   
#   print(b)
# }
# write_csv(data.frame(bstrap1_main), path = "inputs/bootstrap_runs/BootstrapStateEffects_main_us.csv")
# write_csv(data.frame(bstrap1_state), path = "inputs/bootstrap_runs/BootstrapStateEffects_state_us.csv")
# ###########################################################

### US ###

# Load bootstrap runs
bstrap1_main <- read.csv("inputs/bootstrap_runs/BootstrapStateEffects_main_us.csv")[, 1]
bstrap1_state <- read.csv("inputs/bootstrap_runs/BootstrapStateEffects_state_us.csv") %>% as.matrix()

# Load shapefile
usa <- readRDS("inputs/map_boundaries/USA_adm1.rds")

# Load FIPS info
fips <- read.csv("inputs/US_FIPS_Codes.csv", header = TRUE, skip = 1)[, c("State", "FIPS.State")]
fips <- fips[!duplicated(fips), ]
names(fips) <- c("NAME_1", "FIPS")
rownames(fips) <- 1:nrow(fips)
usa@data <- left_join(usa@data, fips, by = "NAME_1") %>% arrange(plotOrder)

# US Regression using feols
mod <- feols(rate_adj ~ tmean * as.factor(state) + prec * as.factor(state) | year, data = data_us, weights = ~popw, cluster = ~fips)

# Pull out coefficients we care about
coef <- coef(mod)
coef <- coef[grep("tmean", names(coef))]
coef <- data.frame(state = as.integer(substr(names(coef), 23, 24)), coef = coef)
coef[1, "state"] <- 1
coef[2:nrow(coef), "coef"] <- coef[2:nrow(coef), "coef"] + coef[1, "coef"] # Everything relative to omitted temp coefficient
rownames(coef) <- 1:nrow(coef)

# Get base rates, population weighted
br <- data_us %>% dplyr::group_by(state) %>% dplyr::summarise(avg = weighted.mean(rate_adj, popw, na.rm = TRUE))

# Calculate % effects = (coef/base rate)*100
coef <- left_join(coef, br, by = "state")

# Get t-stats for significance
st <- sort(unique(data_us$state))
all.equal(rownames(bstrap1_state), as.character(st))
test_stat <- bstrap1_state

for (i in 1:nrow(test_stat)) {
  test_stat[i, ] <- bstrap1_state[i, ] - bstrap1_main
}
pos <- apply(test_stat, 1, function(x) { sum(x > 0) / length(x) })
neg <- apply(test_stat, 1, function(x) { sum(x < 0) / length(x) })

coef$siggy <- (pos > 0.95 | pos < 0.05)

for (s in st) {
  data_us[, paste("temp_x_fe", s, sep = "")] <- as.numeric(data_us$state == s) * data_us$tmean
  data_us[, paste("prec_x_fe", s, sep = "")] <- as.numeric(data_us$state == s) * data_us$prec
}

st <- st[st != 1]

coef$plot <- (coef$coef / coef$avg) * 100
names(coef)[c(1, 5)] <- c("FIPS", "plotValue")
usa@data <- left_join(usa@data, coef[, c("FIPS", "plotValue", "siggy")], by = "FIPS") %>% arrange(plotOrder)

########## PLOT ###########

# Define color palettes
pal.neg <- colorRampPalette(c(brewer.pal(n = 9, name = "Blues")[c(8:5, 1)]))(256)
pal.pos <- colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[c(1, 5:8)]))(256)

breaks.neg <- seq(-4.5, 0, 0.25)
breaks.pos <- seq(0, 7, 0.25)

int.pos.usa <- classIntervals(usa@data$plotValue, style = "fixed", fixedBreaks = breaks.pos)
int.neg.usa <- classIntervals(usa@data$plotValue, style = "fixed", fixedBreaks = breaks.neg)

# Assign colors, handle NA
usa@data$color[usa@data$plotValue > 0 & !is.na(usa@data$plotValue)] <- findColours(int.pos.usa, pal.pos)[usa@data$plotValue > 0 & !is.na(usa@data$plotValue)]
usa@data$color[usa@data$plotValue < 0 & !is.na(usa@data$plotValue)] <- findColours(int.neg.usa, pal.neg)[usa@data$plotValue < 0 & !is.na(usa@data$plotValue)]
usa@data$color[usa@data$plotValue < -4.5 & !is.na(usa@data$plotValue)] <- brewer.pal(n = 9, name = "Blues")[7]
usa@data$color[usa@data$plotValue > 7 & !is.na(usa@data$plotValue)] <- brewer.pal(n = 9, name = "OrRd")[7]

# Handle siggy NA, select significant states
usa_siggy <- usa[usa@data$siggy == TRUE & !is.na(usa@data$siggy), ]

# Generate map
pdf("outputs/raw_figures/Figure5.pdf", width = 6, height = 6)
plot(usa, lwd = 0.025, xlim = c(-125, -66), ylim = c(17, 47), col = usa@data$color, border = NA)
plot(usa_siggy, add = TRUE, lwd = 1.5)
dev.off()