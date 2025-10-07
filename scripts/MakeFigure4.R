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

data_us <- data_us %>%
  mutate(
    date = as.Date(paste(year, month, "01", sep = "-"))
  )

yrs <- rep(2013:2024,each = 12)

mod <- feols(
  rate_adj ~ tmean*factor(date) + prec*factor(date) | fipsmo + stateyear ,
  data = data_us,
  weights = data_us$popw,
  cluster = data_us$fips
)

coef <- summary(mod)$coefficients
coef.se <- summary(mod)$se
vars <- grep("tmean",names(coef))
vcov <- vcov(mod)[vars,vars]

#pull out and proces results  
coef <- coef[vars]
b <- coef[1]
se <- coef.se[1]

for (k in 2:length(vars)) {
  b = c(b,coef[k]+coef[1])
  covv <- vcov[c(1,k),c(1,k)]
  se = c(se,sqrt(sum(covv)))
}

cilo <- b - 1.96*se
cihi <- b + 1.96*se

res <- data.frame(yr=yrs,b,se,cilo,cihi)
avg <- data_us %>% group_by(yr) %>% summarize(avg = weighted.mean(rate_adj,popw))
res <- left_join(res,avg,by="yr")

mod2 <- feols(rate_adj ~ tmean + prec | fipsmo +stateyear , data=data_us, weights=data_us$popw,cluster = data_us$fips)
br <- weighted.mean(data_us$rate_adj,data_us$popw)
eff <- coef(mod2)[1]/br*100

res <- res %>%
  group_by(yr) %>%
  mutate(observation = row_number()) %>%
  ungroup()

res <- res %>% 
  mutate(b = b * 100,
         cilo = cilo * 100,
         cihi = cihi * 100)

p <- ggplot(res) +
  geom_vline(
    aes(xintercept = 0), 
    color = "red", 
    linetype = "dashed", 
    linewidth = 0.7
  ) +
  
  geom_errorbarh(
    aes(xmin = cilo, xmax = cihi, y = factor(observation)),
    height = 0.2, 
    color = "steelblue",
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(x = b, y = factor(observation)),
    size = 1,
    color = "navy"
  ) +
  
  facet_wrap(
    ~ yr, 
    ncol = 3,
    scales = "free_y",
    axes = "all"
  ) +
  
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  
  labs(
    x = "Effect of Monthly Temperature on Police Violence（%change）",
    y = "Month",
    title = ""
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    strip.background = element_rect(fill = "grey80"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(0.8, "lines")
  )

out_dir <- "outputs/raw_figures/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tiff(filename = file.path(out_dir, "Figure4.tiff"),
     width = 8, height = 7.5, units = "in", res = 600)
dev.off()