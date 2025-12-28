library(fect) 
library(dplyr)
library(data.table)
library(zoo)
library(lubridate)
library(fixest)
library(ggplot2)
library(panelView)
library(HonestDiDFEct)

# 1. LOAD DATA
data_use <- read.csv(
  "D:/data.csv",
  header = TRUE
)

str(data_use)

# Filter period up to 6 months after the hurricane
data_use <- data_use %>%
  filter(period < 135) 

data_use <- data_use %>%
  filter(Year > 2012)

# 2. FILTER TREATED AND CONTROL UNITS
data_use <- data_use %>% filter(affected25 %in% c(0,1)) # PLEASE REPLACE WITH affected50 FOR SENSITIVITY ANALYSIS

# 3. BASIC VARIABLES
data_use$id_numeric <- as.integer(as.factor(data_use$DANE))
data_use$period    <- as.integer(data_use$period)
data_use$D <- ifelse(data_use$affected25 == 1 & data_use$treated == 1, 1, 0)   # PLEASE REPLACE WITH affected50 FOR SENSITIVITY ANALYSIS

data_use <- data_use %>% filter(!is.na(cases)) %>% arrange(id_numeric, period)

# 4. CONVERT TO DATA.TABLE
data_use <- as.data.table(data_use)
setkeyv(data_use, c("id_numeric","period"))

# 5. CREATE LAGS OF THE OUTCOME (12)
max_lag <- 12
for(l in 1:max_lag){
  lagname <- paste0("cases_lag",l)
  data_use[, (lagname) := shift(cases, n=l, type="lag"), by=id_numeric]
}

# 6. PRE-TREATMENT SUMMARIES
data_use <- data_use %>%
  group_by(id_numeric) %>%
  mutate(
    mean_cases_pre = lag(cummean(cases)),
    sd_cases_pre   = lag(cummean((cases - mean(cases, na.rm=TRUE))^2)),
    trend_pre      = lag(rollapply(
      cases, width=6,
      FUN=function(x) coef(lm(x ~ seq_along(x)))[2],
      fill=NA, align="right"
    ))
  ) %>%
  ungroup()

# COVARIATES FOR MATCHING
lag_vars     <- paste0("cases_lag", 1:max_lag)

time_varying_vars <- c("t2m","tp","forest")

time_invariant_vars <- c("urban_dim","vectors","log_pop","MPI",
                         "mining","altitude","dist_huracan_km")

pretrend_vars <- c("mean_cases_pre","sd_cases_pre","trend_pre")

all_covariates <- c(lag_vars,
                    pretrend_vars,
                    time_varying_vars,
                    time_invariant_vars)

covs_formula <- as.formula(
  paste("~", paste(all_covariates, collapse = " + "))
)

# MATCHING LAG
lag_num <- 12
cat("Lag used:", lag_num, "\n")

data_use <- as.data.table(data_use)

data_use[, id_numeric := as.integer(id_numeric)]
data_use[, period := as.integer(period)]
data_use[, D := as.integer(D)]
data_use[, cases := as.numeric(cases)]

data_use_df <- as.data.frame(data_use)

# PanelData object
PD <- PanelData(
  panel.data = data_use_df,
  unit.id = "id_numeric",
  time.id = "period",
  treatment = "D",
  outcome = "cases"
)

# MATCHING
PM_results <- PanelMatch(
  PD,
  lag = lag_num,
  refinement.method = "mahalanobis",
  match.missing = FALSE,
  size.match = 5,
  qoi = "att",
  covs.formula = covs_formula,
  use.diagonal.variance.matrix = TRUE
)

# NORMALIZED WEIGHTS + AVERAGE ACROSS TREATED UNITS

matched_weights_list <- list()

for(t_id in names(PM_results$att)){
  
  wts <- attr(PM_results$att[[t_id]], "weights")
  wts <- wts[wts > 0]  # only controls
  
  if(length(wts) > 0){
    
    # normalize within the matched set of each treated unit
    wts_norm <- wts / sum(wts)
    
    matched_weights_list[[t_id]] <- data.frame(
      id_numeric = as.integer(names(wts_norm)),
      weight     = as.numeric(wts_norm)
    )
  }
}

# Average weights across treated units
control_weights_df <- bind_rows(matched_weights_list) %>%
  group_by(id_numeric) %>%
  summarise(weight = mean(weight)) %>%
  ungroup()

# Treated units with weight = 1
treated_ids <- names(PM_results$att) %>% as.integer()
treated_weights_df <- data.frame(id_numeric = treated_ids, weight = 1)

# Combine all weights
all_weights <- bind_rows(control_weights_df, treated_weights_df)

data_matched <- data_use %>% left_join(all_weights, by="id_numeric")
data_matched$weight[is.na(data_matched$weight)] <- 0

cat("Units with weight > 0:", sum(data_matched$weight > 0), "\n")

# fect WITH WEIGHTS

fect_ife_matched <- fect(
  formula = cases ~ D + t2m + tp + forest,
  data = data_matched,
  index = c("id_numeric","period"),
  method = "ife", 
  force = "two-way",
  CV = TRUE,
  r = c(0,5),
  se = TRUE,
  nboots = 2500,
  W = "weight",
  parallel = TRUE,
  criterion = "moment", 
  seed = 123
)

print(fect_ife_matched)
str(fect_ife_matched)
fect_ife_matched$est.avg

# RMSE (Root Mean Square Error) PRE-TREATMENT:
print(rmse_fect <- fect_ife_matched$rmse)

# TABLE 1
tabla_dinamica <- fect_ife_matched$est.att
tabla_dinamica <- as.data.frame(tabla_dinamica)
View(tabla_dinamica)

# FIGURE Observed vs Synthetic Total
cat("Generating Figure 2A: Total Observed vs Synthetic Cases...\n")
# 1. Extract Y (Observed) and Y.ct (Counterfactual) matrices
# --- DATA DIAGNOSTICS ---
cat("--- START DIAGNOSTICS ---\n")
# Check 'period'
cat("Summary of 'period' in data_matched:\n")
print(summary(data_matched$period))
cat("First 20 unique periods:\n")
print(sort(unique(data_matched$period))[1:20])
# Check cases in data_matched for treated units
treated_ids_vec <- unique(data_matched$id_numeric[data_matched$D == 1]) 
cat("Treated IDs (from matching):", paste(head(treated_ids), collapse=","), "...\n")
# Compute total observed cases directly from dataframe
obs_from_df <- data_matched %>%
  filter(id_numeric %in% treated_ids) %>%
  group_by(period) %>%
  summarise(Total_Cases = sum(cases, na.rm = TRUE)) %>%
  arrange(period)
cat("Summary of Total Cases (from DataFrame) for treated units:\n")
print(head(obs_from_df, 20))
cat("...\n")
print(tail(obs_from_df, 20))
# Check for gaps with zero cases
zeros_check <- obs_from_df %>% filter(Total_Cases == 0)
cat("Number of periods with 0 total cases among treated units (DataFrame):", nrow(zeros_check), "\n")
if(nrow(zeros_check) > 0) print(head(zeros_check))
cat("--- END DIAGNOSTICS ---\n")

if (!is.null(fect_ife_matched$Y.dat) && !is.null(fect_ife_matched$Y.ct)) {
  mat_obs <- fect_ife_matched$Y.dat
  mat_syn <- fect_ife_matched$Y.ct
} else if (!is.null(fect_ife_matched$Y) && is.matrix(fect_ife_matched$Y) && !is.null(fect_ife_matched$Y.ct)) {
  mat_obs <- fect_ife_matched$Y
  mat_syn <- fect_ife_matched$Y.ct
} else {
  warning("Y.dat/Y.ct matrices not found in fect object. Attempting reconstruction using Y if available...")
  mat_obs <- NULL
  mat_syn <- NULL
}

# Identify columns corresponding to treated units
col_indices <- NULL

if (!is.null(fect_ife_matched$tr)) {
  cat("Using treatment indices ($tr) from fect object.\n")
  col_indices <- fect_ife_matched$tr
} else {
  # Fallback: try matching by names
  cat("Component $tr not found. Attempting name-based matching...\n")
  unit_names <- colnames(mat_obs)
  if (!is.null(unit_names)) {
    treated_ids_char <- as.character(treated_ids)
    col_indices <- which(unit_names %in% treated_ids_char)
    # If failed, attempt integer matching
    if (length(col_indices) == 0) {
      col_indices <- which(as.integer(unit_names) %in% treated_ids)
    }
  }
}

total_obs <- rowSums(mat_obs[, col_indices, drop = FALSE], na.rm = TRUE)
total_syn <- rowSums(mat_syn[, col_indices, drop = FALSE], na.rm = TRUE)

Observed = total_obs
Synthetic = total_syn
Period <- seq(from = 1, by = 1, to = length(Observed))
index <- seq(from = -127, by = 1, to = 6) 

df_ObsSyn = cbind(Observed, Synthetic, Period, index)

df_ObsSyn <- as.data.frame(df_ObsSyn)

Obs_Syn <- ggplot(df_ObsSyn, aes(x = index)) +
  geom_line(aes(y = Observed, color = "Observed"), size = 1) +
  geom_line(aes(y = Synthetic, color = "Synthetic Control"), size = 0.7) +
  geom_vline(xintercept = df_ObsSyn$index[df_ObsSyn$Period == 128], color = "blue", size = 0.5) + 
  scale_color_manual(name = "Legend", values = c("Observed" = "black", "Synthetic Control" = "red")) +
  labs(
    title = "a",
    y = "Total Cases",
    x = "Periods since the hurricane's onset"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,      
      size = 28,        
      face = "bold"     
    ),
    legend.position = "none",
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

print(Obs_Syn)

# ATT RESULTS
# Extract time-specific treatment effects object
if (!is.null(fect_ife_matched$est.att)) {
  est_att_matrix <- fect_ife_matched$est.att
} else {
  stop("fect_ife_matched$est.att not found. Verify the 'fect' model.")
}

# Convert to data.frame
es_df <- as.data.frame(est_att_matrix)

# Verify conversion
if (!is.data.frame(es_df) || nrow(es_df) == 0) {
  stop("Conversion of 'est.att' to data.frame failed or resulted in an empty dataframe.")
}

# Check original column names
original_names <- colnames(es_df)
print(paste("Original column names:", paste(original_names, collapse = ", ")))

# Rename columns to more manageable names
expected_original_names <- c("ATT", "S.E.", "CI.lower", "CI.upper", "p.value", "count") 
new_names_map <- c("ATT" = "Estimate", "S.E." = "SE", "CI.lower" = "lower95", "CI.upper" = "upper95")
for(old_name in names(new_names_map)){
  if(old_name %in% original_names){
    names(es_df)[names(es_df) == old_name] <- new_names_map[old_name]
  }
}

# Add 'Period' column
es_df$Period <- as.numeric(rownames(es_df))

# Ensure correct data types
numeric_cols <- c("Period", "Estimate", "SE", "lower95", "upper95") 
for(col in numeric_cols){
  if(col %in% names(es_df)){
    es_df[[col]] <- as.numeric(es_df[[col]])
  } else {
    warning(paste("Expected column '", col, "' not found in 'es_df'."))
  }
}

# Confirm required columns exist
required_cols <- c("Period", "Estimate", "SE", "lower95", "upper95")
missing_cols <- setdiff(required_cols, names(es_df))
if(length(missing_cols) > 0) {
  stop("Missing required columns in 'es_df' after processing: ", paste(missing_cols, collapse = ", "))
}

str(es_df)

# Create a monthly date sequence from January 2013 to December 2023
fechas_mensuales <- seq(from = ymd("2013-01-01"), to = ymd("2023-12-31"), by = "1 month")

# Assign calendar year to each row in es_df
es_df <- es_df %>%
  mutate(
    Year_calendar = year(fechas_mensuales[1:nrow(es_df)])
  )

es_df <- es_df %>%
  mutate(Year_calendar = ifelse(is.na(Year_calendar), 2023, Year_calendar))

es_df <- es_df %>%
  mutate(index = seq(from = -127, by = 1, to = 6)) 

# Plot ATT with index on x-axis
p_es <- ggplot(es_df, aes(x = index, y = Estimate)) +
  geom_line(color = "red", size = 1) + 
  geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = es_df$index[es_df$Period == 0], color = "blue", size = 0.5) + 
  labs(
    title = "b",
    y = "Effect on malaria cases",
    x = "Periods since the hurricane's onset"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,      
      size = 28,        
      face = "bold"     
    ),
    legend.position = "none",
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

# Print the plot
print(p_es)

# DIAGNOSTIC TESTS
# Placebo test
out.ife.p <- fect(formula = cases ~ D + t2m + tp + forest,
                  data = data_matched,
                  index = c("id_numeric","period"),
                  method = "ife", 
                  force = "two-way",
                  parallel = TRUE,
                  CV = 0,
                  r = 2,
                  se = TRUE,
                  nboots = 2500,
                  W = "weight",
                  placeboTest = TRUE , placebo.period = c(-6, 0))

plot(out.ife.p, ylab = "Effect of Hurricane Julia on malaria cases", 
     xlab = "Time relative to Hurricane Julia", main = "a", cex.main = 2.5, cex.legend = 1.5, cex.axis = 1.5,
     cex.text = 1.5, cex.lab = 1.5,  stats = c("placebo.p","equiv.p"), show.count = FALSE, xlim = c(-5, 6))

# Tests for (no) pre-trend
plot(fect_ife_matched, type = "equiv", ylim = c(-4,4), ylab = "Effect of Hurricane Julia on malaria cases", 
     xlab = "Time relative to Hurricane Julia", main = "b", cex.main = 2.5, cex.axis = 1.5,
     cex.legend = 1.5,  cex.lab = 1.5, cex.text = 1.5, show.count = FALSE)

# LOO test
out.ife.loo <- fect(cases ~ D + t2m + tp + forest, 
                    data = data_matched, 
                    index = c("id_numeric","period"), 
                    method = "ife", 
                    force = "two-way", 
                    se = TRUE, 
                    parallel = TRUE, 
                    nboots = 2500, 
                    W = "weight",
                    loo = TRUE)

plot(out.ife.loo, type = "equiv", ylim = c(-5,6), loo = TRUE, ylab = "Effect of Hurricane Julia on malaria cases",
     xlab = "Time relative to Hurricane Julia", cex.lab = 1.5, cex.axis = 1.5,
     cex.legend = 1.5, main = "c", cex.main = 2.5, cex.text = 1.5,
     show.count = FALSE)

# Conditional Average Treatment Effect on the Treated (CATT)
plot(fect_ife_matched, type = "hte", covariate = "t2m", ylab = "Effect of Hurricane Julia on malaria cases",
     xlab = "Temperature (°C)", main = "a", cex.lab = 1.5, cex.axis = 1.5,
     cex.legend = 1.5, cex.main = 2.5, cex.text = 1.5)
plot(fect_ife_matched, type = "hte", covariate = "tp", ylab = "Effect of Hurricane Julia on malaria cases",
     xlab = "Rainfall (mm)", main = "b", cex.lab = 1.5, cex.axis = 1.5,
     cex.legend = 1.5, cex.main = 2.5, cex.text = 1.5)
plot(fect_ife_matched, type = "hte", covariate = "forest",
     ylab = "Effect of Hurricane Julia on malaria cases",
     xlab = "Forest coverage (%)", main = "c", cex.lab = 1.5, cex.axis = 1.5,
     cex.legend = 1.5, cex.main = 2.5, cex.text = 1.5)

# SENSITIVITY ANALYSIS: RELATIVE MAGNITUDE (RM) RESTRICTION
T.post <- 6 # Number of post-treatment periods based on original analysis
post_periods_vec <- 1:T.post

# Parameters for Relative Magnitude (RM) restriction
Mbar_vec_avg_rm <- seq(0, 1, by = 0.1)    
Mbar_vec_period_rm <- c(0, 0.1)          

# Parameters for Smoothness restriction
M_vec_avg_smooth <- seq(0, 0.25, by = 0.05) 
M_vec_period_smooth <- c(0, 0.1)           

out.fect.placebo <- fect_sens(
  fect.out      = out.ife.p,
  post.periods  = post_periods_vec,
  Mbarvec       = Mbar_vec_avg_rm,
  periodMbarvec = Mbar_vec_period_rm,
  Mvec          = M_vec_avg_smooth,
  periodMvec    = M_vec_period_smooth,
  parallel      = TRUE 
)

plot(out.fect.placebo,
     type = "sens",
     restrict = "rm",
     main = "Relative Magnitude Restriction")

plot(out.fect.placebo,
     type = "sens_es",
     restrict = "rm",
     main = "ATTs with Robust Confidence Sets (RM)",
     ylab = "Coefficients and 95% CI",
     xlim = c(-12,6), 
     ylim = c(-8,6), 
     show.count = FALSE,
     sens.colors = c("blue", "red"))

# SENSITIVITY ANALYSIS: SMOOTHNESS RESTRICTION
plot(out.fect.placebo,
     type = "sens",
     restrict = "sm",
     main = "Smoothness Restriction")

plot(out.fect.placebo,
     type = "sens_es",
     restrict = "sm",
     main = "ATTs with Robust Confidence Sets (Smoothness)",
     ylab = "Coefficients and 95% CI",
     xlim = c(-12,6),
     ylim = c(-15,5),
     show.count = FALSE,
     sens.colors = c("blue", "red"))

