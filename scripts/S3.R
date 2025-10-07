file_path <- "inputs/data/All/Mapping Police Violence_new"
dta <- fread(file_path)

setorder(dta, fips, stateyear)

dta[, tmean := scale(tmean)]
dta[, prec := scale(prec)]

dta[, tmean_lag0 := tmean]
dta[, prec_lag0 := prec]
for (i in 1:2) {
  dta[, paste0("tmean_lag", i) := shift(tmean, i, type = "lag"), by = fips]
  dta[, paste0("tmean_lead", i) := shift(tmean, i, type = "lead"), by = fips]
  dta[, paste0("prec_lag", i) := shift(prec, i, type = "lag"), by = fips]
  dta[, paste0("prec_lead", i) := shift(prec, i, type = "lead"), by = fips]
}

dta[, fipsmo := paste0(fips, "-", sprintf("%02d", month))]

vars_to_check <- c("rate_adj", "tmean_lag0", paste0("tmean_lag", 1:2), paste0("tmean_lead", 1:2),
                   "prec_lag0", paste0("prec_lag", 1:2), paste0("prec_lead", 1:2),
                   "popw", "stateyear", "fipsmo", "fips")
dta <- dta[complete.cases(dta[, ..vars_to_check]), ]

runreg_mod <- function(lag = 2, lead = 2, y = "rate_adj", temp = "tmean", prec = "prec",
                       data, weights = data$popw, fe = "fipsmo+stateyear", cl = "fips") {
  
  if (lead > 0) {
    temp_vars <- c(paste0(temp, "_lead", lead:1), paste0(temp, "_lag", 0:lag))
    prec_vars <- c(paste0(prec, "_lead", lead:1), paste0(prec, "_lag", 0:lag))
  } else {
    temp_vars <- paste0(temp, "_lag", 0:lag)
    prec_vars <- paste0(prec, "_lag", 0:lag)
  }
  
  all_vars <- c(temp_vars, prec_vars)
  lgs <- paste(all_vars, collapse = " + ")
  
  fmla <- as.formula(paste0("rate_adj ~ ", lgs, " | stateyear | 0 | fips"))
  mod <- felm(fmla, data = data, weights = weights)
  
  coef_mod <- summary(mod)$coefficients
  mn <- weighted.mean(data[[y]], weights, na.rm = TRUE)
  
  temp_coef <- coef_mod[temp_vars, , drop = FALSE]
  
  temp_group <- if (lead > 0) {
    c(paste0("t+", lead:1), "t", paste0("t-", 1:lag))
  } else {
    c("t", paste0("t-", 1:lag))
  }
  
  temp_r1 <- data.frame(
    est   = temp_coef[, 1],
    se    = temp_coef[, 2],
    mean  = mn,
    var   = "Temperature",
    group = temp_group,
    N     = mod$N
  )
  
  lag_indices <- if (lead > 0) {
    (lead + 1):(lead + lag + 1)
  } else {
    1:(lag + 1)
  }
  
  temp_r2 <- data.frame(
    est   = sum(temp_coef[lag_indices, 1]),
    se    = sqrt(sum(vcov(mod)[temp_vars[lag_indices], temp_vars[lag_indices]])),
    mean  = mn,
    var   = "Temperature",
    group = "combined",
    N     = mod$N
  )
  
  out <- rbind(temp_r1, temp_r2)
  
  out$ci_lower <- out$est - 1.96 * out$se  
  out$ci_upper <- out$est + 1.96 * out$se  
  
  out <- as.data.table(out)
  return(out)
}

reg_result <- runreg_mod(lag = 2, lead = 2, data = dta)

reg_result$group <- factor(reg_result$group, levels = c("combined", paste0("t-", 2:1), "t", paste0("t+", 1:2)))

p <- ggplot(reg_result, aes(x = group, color = var)) +
  geom_point(aes(y = est * 100), position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = (est - se) * 100, ymax = (est + se) * 100),
                width = 0.2,
                position = position_dodge(width = 0.3)) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  theme_bw() +
  labs(
    x = "Time Period",
    y = "% temperature effect on death rate of police violence"
  ) +
  scale_color_manual(values = c("Temperature" = "red")) +
  theme(
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = "white"),
    axis.text.y           = element_text(angle = 45, hjust = 1,size = 8)
    
  )

pdf_file <- "outputs/raw_figures/S3.pdf"
ggsave(pdf_file, plot = p, width = 5, height = 5, dpi = 600)
cat("ok", pdf_file, "\n")