source("scripts/functions.R")
setwd("scripts")
data_us <- read.csv("inputs/data/All/Mapping Police Violence_new")
data_us <- data_us %>% dplyr::filter(year >= 2022 & year <= 2024) #The value can be adjusted within the specified range to match the desired year
data_us <- data_us[-which(data_us$state == 2),] 

usa <- readRDS("inputs/map_boundaries/USA_adm1.rds")
library(sf)
usa <- st_as_sf(usa)

area <- readxl::read_xlsx("inputs/NOAA_Climate_Regions_2024.xlsx")
fips <- read.csv("US_FIPS_Codes.csv",skip = 1)
head(fips)
fips <- fips[,c(1,3)]
fips <- unique(fips)
dim(fips)

colnames(fips) <- c("State", "FIPS")
area <- merge(area, fips, by = "State")

area <- merge(area,fips,by = "State")
colnames(area)[3] <- "state"
dim(area)

head(area)
head(data_us)
data_us <- merge(data_us,area,by = "state")
dim(data_us) 

fips_from_geom <- data.frame(State = usa$NAME_1,fips_state = usa$ID_1)
fips_from_geom <- fips_from_geom[-which(fips_from_geom$State == "District of Columbia"),]
identical(fips_from_geom$State,unique(data_us$State))

data_us <- merge(data_us,fips_from_geom,by = "State")
dim(data_us) 
head(data_us)
data_us <- data_us %>% dplyr::select(-c("state"))

###########################################################
# The below code has been commented out because it
# takes several hours to run. The stored data
# from these bootstrap runs will be loaded below.
# Uncomment this section to run bootstraps and replace
# file `BootstrapStateEffects_main_us.csv` and `BootstrapStateEffects_state_us.csv` with updated data.
###########################################################
B <- 1000
n1 <- nrow(data_us)

data_us$bstrap_id <- 1:nrow(data_us)

bstrap1_main <- rep(NA, B); names(bstrap1_main) <- paste("b", 1:B, sep = "")
bstrap1_state <- matrix(nrow = length(unique(data_us$Climate_Region)), ncol = B) # 9
rownames(bstrap1_state) <- sort(unique(data_us$Climate_Region))
colnames(bstrap1_state) <- paste("b", 1:B, sep = "")

us_states <- data.frame(state = sort(unique(data_us$Climate_Region)))

# Bootstrap
# for (b in 1:B) {
#   data_us.sample <- left_join(data.frame(bstrap_id = sample(x = data_us$bstrap_id, size = n1, replace = TRUE)), data_us, "bstrap_id")
#   while(length(unique(data_us.sample$Climate_Region)) < 9) {
#     data_us.sample <- left_join(data.frame(bstrap_id = sample(x = data_us$bstrap_id, size = n1, replace = TRUE)), data_us, "bstrap_id")
#   }
# 
#   # US model using feols
#   mod1_main <- feols(rate_adj ~ tmean + prec | Climate_Region + stateyear, data = data_us.sample, weights = ~popw)
#   mod1_state <- feols(rate_adj ~ tmean * as.factor(Climate_Region) + prec * as.factor(Climate_Region) | fipsmo + stateyear, data = data_us.sample, weights = ~popw)
# 
#   bstrap1_main[b] <- coef(mod1_main)["tmean"]
#   coef <- coef(mod1_state)[grep("tmean", names(coef(mod1_state)))]
#   coef[2:length(coef)] <- coef[1] + coef[2:length(coef)]
#   bstrap1_state[, b] <- as.numeric(coef)
# 
# }
# 
# write_csv(data.frame(bstrap1_main), path = "inputs/bootstrap_runs/Bootstrapclimate_regionEffects_main_us_all.csv")
# write_csv(data.frame(bstrap1_state), path = "inputs/bootstrap_runs/Bootstrapclimate_regionEffects_state_us_all.csv")
# ###########################################################

### US ###

# Load bootstrap runs
bstrap1_main <- read.csv("inputs/bootstrap_runs/Bootstrapclimate_regionEffects_main_us_all.csv")[, 1]
bstrap1_state <- read.csv("inputs/bootstrap_runs/Bootstrapclimate_regionEffects_state_us_all.csv") %>% as.matrix()


# Load shapefile
usa <- readRDS("inputs/map_boundaries/USA_adm1.rds")
library(sf)
usa <- st_as_sf(usa)
# Load FIPS info
usa <- usa[-which(usa$OBJECTID == 10),]



names(area)[1] <- c("NAME_1")
rownames(area) <- 1:nrow(area)
usa <- left_join(usa, area, by = c("NAME_1")) 

usa_border <- usa %>%  
  group_by(Climate_Region) %>% summarise(geometry = st_union(geometry)) %>% 
  st_simplify(dTolerance = 10) 

plot(usa_border)

usa_map_region <- usa_border %>%
  group_by(Climate_Region) %>% 
  mutate(
    centroid = st_centroid(geometry),
    lon = st_coordinates(centroid)[, "X"],
    lat = st_coordinates(centroid)[, "Y"]
  )


# US Regression using felm
mod <- feols(rate_adj ~ tmean * as.factor(Climate_Region) + prec * as.factor(Climate_Region) | stateyear, data = data_us, weights = ~popw, cluster = ~Climate_Region)

# Pull out coefficients we care about
coef <- coef(mod)
coef <- coef[grep("tmean", names(coef))]
coef <- data.frame(Climate_Region = substr(names(coef), 32, nchar(names(coef))), coef = coef)
coef[1, "Climate_Region"] <- "Northeast"  #Northeast
coef[2:nrow(coef), "coef"] <- coef[2:nrow(coef), "coef"] + coef[1, "coef"] # Everything relative to omitted temp coefficient
rownames(coef) <- 1:nrow(coef)

# Get base rates, population weighted
br <- data_us %>% dplyr::group_by(Climate_Region) %>% dplyr::summarise(avg = weighted.mean(rate_adj, popw, na.rm = TRUE))

# Calculate % effects = (coef/base rate)*100
coef <- left_join(coef, br, by = "Climate_Region")

# Get t-stats for significance
st <- sort(unique(data_us$Climate_Region))
all.equal(rownames(bstrap1_state), as.character(st))
test_stat <- bstrap1_state

for (i in 1:nrow(test_stat)) {
  test_stat[i, ] <- bstrap1_state[i, ] - bstrap1_main
}
pos <- apply(test_stat, 1, function(x) { sum(x > 0) / length(x) })
neg <- apply(test_stat, 1, function(x) { sum(x < 0) / length(x) })

coef$siggy <- (pos > 0.95 | pos < 0.05)

for (s in st) 
{
  data_us[, paste("temp_x_fe", s, sep = "")] <- as.numeric(data_us$Climate_Region == s) * data_us$tmean
  data_us[, paste("prec_x_fe", s, sep = "")] <- as.numeric(data_us$Climate_Region == s) * data_us$prec
}

#st <- st[st != 1]

coef$plot <- (coef$coef / coef$avg) * 100


names(coef)[c(1, 5)] <- c("Climate_Region", "plotValue")
usa <- left_join(usa_border, coef[, c("Climate_Region", "plotValue", "siggy")], by = "Climate_Region") 

########## PLOT ###########
library(RColorBrewer)
# Define color palettes
pal.neg <- colorRampPalette(c(brewer.pal(n = 9, name = "Blues")[c(8:5, 1)]))(256)
pal.pos <- colorRampPalette(c(brewer.pal(n = 9, name = "Reds")[c(1, 5:8)]))(256)

breaks.neg <- seq(-0.2, 0, 0.1)
breaks.pos <- seq(0, 1.2, 0.1)

library(classInt)
int.pos.usa <- classIntervals(usa$plotValue, style = "fixed", fixedBreaks = breaks.pos)
int.neg.usa <- classIntervals(usa$plotValue, style = "fixed", fixedBreaks = breaks.neg)

# Assign colors, handle NA
usa$color[usa$plotValue > 0 & !is.na(usa$plotValue)] <- findColours(int.pos.usa, pal.pos)[usa$plotValue > 0 & !is.na(usa$plotValue)]
usa$color[usa$plotValue < 0 & !is.na(usa$plotValue)] <- findColours(int.neg.usa, pal.neg)[usa$plotValue < 0 & !is.na(usa$plotValue)]
usa$color[usa$plotValue < -0.1 & !is.na(usa$plotValue)] <- brewer.pal(n = 9, name = "Blues")[7]
usa$color[usa$plotValue > 0.8 & !is.na(usa$plotValue)] <- brewer.pal(n = 9, name = "OrRd")[7] 

# Handle siggy NA, select significant states
usa_siggy <- usa[usa$siggy == TRUE & !is.na(usa$siggy), ]

# Generate map
pdf("outputs/raw_figures/Figure6.pdf", width = 6, height = 6)
plot(usa$geometry, lwd = 0.025, xlim = c(-125, -66), ylim = c(17, 47), col = usa$color, border = NA)
plot(usa_siggy$geometry, add = TRUE, lwd = 1.5)
text(x = usa_map_region$lon,
     y = usa_map_region$lat,
     labels = usa_map_region$Climate_Region,  
     adj = 0.5,
     cex = 0.52,            
     col = "black",        
     font = 1.5)            
dev.off()