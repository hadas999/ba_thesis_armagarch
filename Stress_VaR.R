# 0. Load libraries ---------------------------------------------------------
library(quantmod)
library(rugarch)
library(PerformanceAnalytics)
library(FinTS)
library(xts)
library(zoo)
library(moments)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(gridExtra)
library(sandwich)
library(aod)
library(lmtest)
library(car)
# 1. Data & Stress Window --------------------------------------------------
ivv <- getSymbols("IVV", src="yahoo", from="2000-01-01", to="2025-05-28", auto.assign=FALSE)
r   <- na.omit(diff(log(Ad(ivv))))
stress <- r["2007-01-01/2009-03-31"]

# 2. Specify base ARMA(1,1)-GARCH(1,1) ------------------------------------
base_spec <- list(
  mean.model     = list(armaOrder=c(1,1), include.mean=TRUE),
  variance.model = list(model="sGARCH", garchOrder=c(1,1))
)

# 3. Fit on Stress Period ---------------------------------------------------
spec_norm_stress <- ugarchspec(mean.model=base_spec$mean.model,
                               variance.model=base_spec$variance.model,
                               distribution.model="norm")
spec_t_stress    <- ugarchspec(mean.model=base_spec$mean.model,
                               variance.model=base_spec$variance.model,
                               distribution.model="std")
fit_norm <- ugarchfit(spec_norm_stress, data=stress, solver="hybrid")
fit_t    <- ugarchfit(spec_t_stress,    data=stress, solver="hybrid")

# 4. Extract & Lock Parameters ---------------------------------------------
pars_norm <- coef(fit_norm)
pars_t    <- coef(fit_t)
nu        <- pars_t["shape"]

spec_norm_full <- ugarchspec(mean.model=base_spec$mean.model,
                             variance.model=base_spec$variance.model,
                             distribution.model="norm",
                             fixed.pars=as.list(pars_norm))
spec_t_full    <- ugarchspec(mean.model=base_spec$mean.model,
                             variance.model=base_spec$variance.model,
                             distribution.model="std",
                             fixed.pars=as.list(pars_t))

# 5. Filter Full Series & Compute Forecasts --------------------------------
f_norm <- ugarchfilter(spec_norm_full, data=r)
f_t    <- ugarchfilter(spec_t_full,    data=r)
mu_n  <- fitted(f_norm);    sd_n  <- sigma(f_norm)
mu_t  <- fitted(f_t);       sd_t  <- sigma(f_t)

# 6. Compute VaR & ES for 95% & 99% -----------------------------------------
alphas <- c(0.05, 0.01)
df <- data.frame(Date=index(r), Return=as.numeric(r))
for(alpha in alphas) {
  qn     <- qnorm(alpha)
  df[[paste0("VaR",alpha*100,"_norm")]] <- mu_n + sd_n * qn
  df[[paste0("ES", alpha*100,"_norm")]] <- mu_n - sd_n*(dnorm(qn)/alpha)
  qtval  <- qt(alpha, df=nu)
  df[[paste0("VaR",alpha*100,"_t")]] <- mu_t + sd_t * qtval
  df[[paste0("ES", alpha*100,"_t")]] <- mu_t - sd_t*(dt(qtval,df=nu)*(nu+qtval^2)/((nu-1)*alpha))
}

# mark violations
for(alpha in alphas) {
  pct <- alpha*100
  df[[paste0("vio",pct,"_norm")]] <- ifelse(df$Return < df[[paste0("VaR",pct,"_norm")]], df$Return, NA)
  df[[paste0("vio",pct,"_t")]]    <- ifelse(df$Return < df[[paste0("VaR",pct,"_t")]],    df$Return, NA)
}
#mark violations for ES
for(alpha in alphas) {
  pct <- alpha*100
  df[[paste0("vioE",pct,"_norm")]] <- ifelse(df$Return < df[[paste0("ES",pct,"_norm")]], df$Return, NA)
  df[[paste0("vioE",pct,"_t")]]    <- ifelse(df$Return < df[[paste0("ES",pct,"_t")]],    df$Return, NA)
}

# 7. Plot VaR & ES Separately ------------------------------------------------
plot_vars <- function(alpha, dist) {
  # alpha = 0.05 or 0.01
  pct      <- alpha * 100
  var_col  <- paste0("VaR", pct, "_", dist)
  es_col   <- paste0("ES",  pct, "_", dist)
  vio_col  <- paste0("vio", pct, "_", dist)
  dist_lbl <- if(dist == "norm") "Normal" else "Student’s t"
  
  ggplot(df, aes(x = Date)) +
    # actual returns
    geom_line(aes(y = Return), linewidth = 0.3, colour = "grey40") +
    # VaR
    geom_line(aes(y = .data[[var_col]]),
              color    = if(dist == "norm") "red" else "blue",
              linewidth = 0.5) +
    # ES
    #geom_line(aes(y = .data[[es_col]]),
    #         color    = if(dist == "norm") "darkred" else "darkblue",
    #        linetype = "dashed", linewidth = 0.8) +
    # Violation points
    geom_point(aes(y = .data[[vio_col]]),
               color = "black", shape = 17, size = 2.5, na.rm = TRUE) +
    labs(
      title = sprintf("VaR at %.0f%% (%s)", 100-pct, dist_lbl),
      y     = "Log Return / Risk Measure"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}
# plot function for ES with violations
plot_es <- function(alpha, dist) {
  # alpha = 0.05 or 0.01
  pct      <- alpha * 100
  es_col   <- paste0("ES",  pct, "_", dist)
  vio_col  <- paste0("vioE", pct, "_", dist)
  dist_lbl <- if(dist == "norm") "Normal" else "Student’s t"
  
  ggplot(df, aes(x = Date)) +
    # actual returns
    geom_line(aes(y = Return), linewidth = 0.3, colour = "grey40") +
    # ES
    geom_line(aes(y = .data[[es_col]]),
              color    = if(dist == "norm") "red" else "blue",
              linewidth = 0.5) +
    # Violation points
    geom_point(aes(y = .data[[vio_col]]),
               color = "black", shape = 17, size = 2.5, na.rm = TRUE) +
    labs(
      title = sprintf("Expected Shortfall at %.0f%% (%s)", 100-pct, dist_lbl),
      y     = "Log Return / Risk Measure"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}



# Now call it like this:
print(plot_vars(0.05, "norm"))  # 95% Normal
print(plot_vars(0.05, "t"))     # 95% Student’s t
print(plot_vars(0.01, "norm"))  # 99% Normal
print(plot_vars(0.01, "t"))     # 99% Student’s t
# Plot ES
print(plot_es(0.05, "norm"))  # 95% Normal ES
print(plot_es(0.05, "t"))     # 95% Student’s t ES
print(plot_es(0.01, "norm"))  # 99% Normal ES
print(plot_es(0.01, "t"))     # 99% Student’s t ES

#number of violations
num_violations <- function(alpha, dist) {
  pct <- alpha * 100
  vio_col <- paste0("vio", pct, "_", dist)
  sum(!is.na(df[[vio_col]]))
}
# number of expected violations
expected_violations <- function(alpha) {
  pct <- alpha * 100
  nrow(df) * alpha
}
#print number of violations
cat("Number of Violations (95% Normal):", num_violations(0.05, "norm"), "\n")
cat("Number of Violations (95% t):", num_violations(0.05, "t"), "\n")
cat("Number of Violations (99% Normal):", num_violations(0.01, "norm"), "\n")
cat("Number of Violations (99% t):", num_violations(0.01, "t"), "\n")
cat("Expected Violations (95%):", expected_violations(0.05), "\n")
cat("Expected Violations (99%):", expected_violations(0.01), "\n")

#number of violations ES
num_violations_es <- function(alpha, dist) {
  pct <- alpha * 100
  vio_col <- paste0("vioE", pct, "_", dist)
  sum(!is.na(df[[vio_col]]))
}
# number of expected violations ES
expected_violations_es <- function(alpha) {
  pct <- alpha * 100
  nrow(df) * alpha
}
#print number of violations ES
cat("Number of Violations (99% Normal):", num_violations_es(0.01, "norm"), "\n")
cat("Number of Violations (99% t):", num_violations_es(0.01, "t"), "\n")
cat("Expected Violations (99%):", expected_violations_es(0.01), "\n")


# 8. Residual Diagnostics ---------------------------------------------------
# Residuals
res_n <- residuals(fit_norm, standardize = TRUE)
res_t <- residuals(fit_t,    standardize = TRUE)

# Histograms + density
p_hist_n <- ggplot(data.frame(res=res_n), aes(res)) +
  geom_histogram(aes(y=..density..), bins=50, fill="grey80", color="black") +
  stat_function(fun=dnorm, args=list(mean=0,sd=1), color="red") +
  labs(title="Std Residuals Histogram (Normal)")

p_hist_t <- ggplot(data.frame(res=res_t), aes(res)) +
  geom_histogram(aes(y=..density..), bins=50, fill="grey80", color="black") +
  stat_function(fun=dt, args=list(df=nu), color="blue") +
  labs(title=sprintf("Std Residuals Histogram (t, ν=%.1f)", nu))

# QQ-Plots
p_qq_n <- ggplot(data.frame(res=res_n), aes(sample=res)) +
  stat_qq() + stat_qq_line(color="red") +
  labs(title="QQ Plot: Normal Residuals")

p_qq_t <- ggplot(data.frame(res=res_t), aes(sample=res)) +
  stat_qq(distribution=qt, dparams=list(df=nu)) +
  stat_qq_line(color="blue") +
  labs(title="QQ Plot: t Residuals")

grid.arrange(p_hist_n, p_hist_t, p_qq_n, p_qq_t, ncol=2)

mse_n <- mean(res_n^2)
    
mse_t <- mean(res_t^2)
# Print MSE
print(sprintf("MSE (Normal): %.4f", mse_n))
print(sprintf("MSE (t): %.4f", mse_t))


# 9. Backtesting (Kupiec & Christoffersen) -----------------------------------
  
vars      <- c("VaR5_norm","VaR5_t","VaR1_norm","VaR1_t")
alphas    <- c(0.05,       0.05,     0.01,       0.01)
res_list  <- vector("list", length(vars))


for(i in seq_along(vars)) {
  var_col      <- vars[i]
  alpha_target <- alphas[i]
  
  df_tmp <- df %>%
    arrange(Date) %>%
    mutate(
      I     = as.integer(Return < .data[[var_col]]),
      I_lag = lag(I)
    ) %>%
    filter(!is.na(I_lag))   # only drop the first NA lag, keep I=0 and I=1
  
  # 1) Kupiec UC
  mod_uc <- lm(I ~ 1, data = df_tmp)
  uc_ct  <- coeftest(mod_uc, vcov = vcovHC(mod_uc, "HC1"))
  uc_est <- uc_ct["(Intercept)", "Estimate"]
  uc_se  <- uc_ct["(Intercept)", "Std. Error"]
  uc_t   <- (uc_est - alpha_target) / uc_se
  uc_p   <- 2 * pt(-abs(uc_t), df = df.residual(mod_uc))
  
  # 2) Christoffersen IND
  mod_ind  <- lm(I ~ 1 + I_lag, data = df_tmp)
  ind_ct   <- coeftest(mod_ind, vcov = vcovHC(mod_ind, "HC1"))
  beta_est <- ind_ct["I_lag","Estimate"]
  beta_se  <- ind_ct["I_lag","Std. Error"]
  beta_t   <- beta_est / beta_se
  beta_p   <- ind_ct["I_lag","Pr(>|t|)"]
  #-----------------------------------------------
  mod <- lm(I ~ I_lag, data = df_tmp)
  
  # 3) Joint Wald test via car::linearHypothesis
  # ---------------------------------------------
  # install.packages("car"); install.packages("sandwich")
  #library(car)
  #library(sandwich)
  
  # Define the two linear hypotheses:
  #   (Intercept) = alpha   AND   I_lag = 0
  L <- rbind(
    "(Intercept) = a" = c(1, 0),  # 1*beta0 + 0*beta1 = alpha
    "I_lag = 0"       = c(0, 1)   # 0*beta0 + 1*beta1 = 0
  )
  # The RHS vector must match in order:
  rhs <- c(alpha = alpha_target, beta = 0)
  
  # Choose a Newey–West lag (e.g. the default Andrews rule of thumb)
  nw_lag <- floor(4 * (nrow(df) / 100)^(2/9))
  
  # Compute the Newey–West HAC covariance matrix
  V_hac <- NeweyWest(mod,
                     lag       = nw_lag,
                     prewhite  = FALSE,
                     adjust    = TRUE)
  
  # Perform the robust Wald test
  lh <- linearHypothesis(
    model = mod,
    hypothesis.matrix = L,
    rhs = rhs,
    # use Newey-West covariance matrix
    vcov = V_hac,
    test = "Chisq"
  )
  
  res_list[[i]] <- data.frame(
    VaR_Series = var_col,
    Alpha      = alpha_target,
    Test       = c("Kupiec_UC", "Christoffersen_IND", "Conditional_CC"),
    Statistic  = c(uc_t, beta_t, lh$Chisq[2]),
    p_value    = c(uc_p, beta_p, lh$`Pr(>Chisq)`[2]),
    stringsAsFactors = FALSE
  )
  
}

results_df <- do.call(rbind, res_list)
print(results_df)
