#───────────────────────────────────────────────────────────────────────────────
# 0. Install & load packages + define backtest functions
#───────────────────────────────────────────────────────────────────────────────
library(quantmod); library(rugarch); library(tseries)
library(PerformanceAnalytics); library(FinTS)
library(xts); library(zoo); library(moments)
library(ggplot2); library(gridExtra)
library(lubridate); library(reshape2)
library(dplyr); library(tidyr)
library(aod); library(lmtest) 
library(sandwich); library(car)
 #───────────────────────────────────────────────────────────────────────────────
 # 1. Data Download & Log-Return Preprocessing
 #───────────────────────────────────────────────────────────────────────────────
 ivv_xts <- getSymbols("IVV", src="yahoo", from="2000-01-01", to="2025-05-28",
                       auto.assign=FALSE)
 prices  <- Ad(ivv_xts)
 r       <- na.omit(diff(log(prices)))
 
 # Quick diagnostic
 cat("ADF p-value:", adf.test(r)$p.value, "\n")
 cat("Skewness:", skewness(r), "Kurtosis:", kurtosis(r), "\n")
 print(jarque.bera.test(r))
 
 # Plot the log-returns
 ggplot(data = data.frame(Date = index(r), Return = as.numeric(r)),
        aes(x = Date, y = Return)) +
   geom_line(color = "grey40", size = 0.5) +
   labs(title = "Log-Returns of IVV ETF",
        x = "Date", y = "Log-Return") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5))
 #───────────────────────────────────────────────────────────────────────────────
 # 2. Rolling-Window Forecast Setup
 #───────────────────────────────────────────────────────────────────────────────
 W     <- 1000
 N     <- length(r)
 steps <- N - W
 alphas <- c(0.05, 0.01)   # 5% and 1%
 
 # Preallocate arrays: rows=steps, cols=distributions×alphas
 VaR_norm <- ES_norm <- array(NA, dim=c(steps,2),
                              dimnames=list(NULL, paste0(c("5","1"),"%")))
 VaR_t    <- ES_t    <- VaR_norm
 actuals  <- rep(NA, steps)
 mu_hat   <- rep(NA, steps)  # for MSE
 #───────────────────────────────────────────────────────────────────────────────
 # 3. Define ARMA(1,1)-GARCH(1,1) specs
 #───────────────────────────────────────────────────────────────────────────────
 spec_norm <- ugarchspec(mean.model=list(armaOrder=c(1,1),include.mean=TRUE),
                         variance.model=list(model="sGARCH",garchOrder=c(1,1)),
                         distribution.model="norm")
 spec_t    <- ugarchspec(mean.model=list(armaOrder=c(1,1),include.mean=TRUE),
                         variance.model=list(model="sGARCH",garchOrder=c(1,1)),
                         distribution.model="std")
 #───────────────────────────────────────────────────────────────────────────────
 # 4. Rolling Forecast Loop
 #───────────────────────────────────────────────────────────────────────────────
 for(i in 1:steps) {
   window_data <- r[i:(i+W-1)]
   next_obs    <- as.numeric(r[i+W])
   
   fit_n <- tryCatch(ugarchfit(spec_norm, data=window_data,
                               solver="hybrid", silent=TRUE),
                     error=function(e)NULL)
   fit_s <- tryCatch(ugarchfit(spec_t, data=window_data,
                               solver="hybrid", silent=TRUE),
                     error=function(e)NULL)
   if(is.null(fit_n)||is.null(fit_s)) next
   
   fc_n <- ugarchforecast(fit_n, n.ahead=1)
   fc_s <- ugarchforecast(fit_s, n.ahead=1)
   
   mu_n <- as.numeric(fitted(fc_n));  sd_n <- as.numeric(sigma(fc_n))
   mu_s <- as.numeric(fitted(fc_s));  sd_s <- as.numeric(sigma(fc_s))
   nu   <- coef(fit_s)["shape"]
   
   mu_hat[i] <- mu_s  # mean forecast
   
   for(j in 1:2) {
     α <- alphas[j]
     # Normal
     qn <- qnorm(α)
     VaR_norm[i,j] <- mu_n + sd_n*qn
     ES_norm[i,j]  <- mu_n - sd_n*(dnorm(qn)/α)
     # Student t
     qtj <- qt(α, df=nu)
     VaR_t[i,j]    <- mu_s + sd_s*qtj
     ES_t[i,j]     <- mu_s - sd_s*(dt(qtj,df=nu)*(nu+qtj^2)/((nu-1)*α))
   }
   actuals[i] <- next_obs
 }
 
 dates <- index(r)[(W+1):N]
 df <- data.frame(Date=dates, Return=actuals,
                  VaR_norm,ES_norm,VaR_t,ES_t)
 colnames(df)[3:6] <- c("VaR5_norm","VaR1_norm","ES5_norm","ES1_norm")
 colnames(df)[7:10]<- c("VaR5_t","VaR1_t","ES5_t","ES1_t")
 
 # 1) Build the data.frame
 dates <- index(r)[(W+1):N]
 df <- data.frame(
   Date     = dates,
   Return   = actuals,
   VaR5_norm = VaR_norm[,1],
   VaR1_norm = VaR_norm[,2],
   ES5_norm  = ES_norm[,1],
   ES1_norm  = ES_norm[,2],
   VaR5_t    = VaR_t[,1],
   VaR1_t    = VaR_t[,2],
   ES5_t     = ES_t[,1],
   ES1_t     = ES_t[,2]
 )
 
 # 2) Rename to consistent names (if needed)
 # colnames(df)[3:6]  <- c("VaR5_norm","VaR1_norm","ES5_norm","ES1_norm")
 # colnames(df)[7:10] <- c("VaR5_t","VaR1_t","ES5_t","ES1_t")
 
 # 3) Add violation flags
 df$vio5_norm <- ifelse(df$Return < df$VaR5_norm, df$Return, NA)
 df$vio1_norm <- ifelse(df$Return < df$VaR1_norm, df$Return, NA)
 df$vio5_t    <- ifelse(df$Return < df$VaR5_t,    df$Return, NA)
 df$vio1_t    <- ifelse(df$Return < df$VaR1_t,    df$Return, NA)
 df$es5_norm  <- ifelse(df$Return < df$ES5_norm,  df$Return, NA)
 df$es1_norm  <- ifelse(df$Return < df$ES1_norm,  df$Return, NA)
 df$es5_t     <- ifelse(df$Return < df$ES5_t,     df$Return, NA)
 df$es1_t     <- ifelse(df$Return < df$ES1_t,     df$Return, NA)
 
 #───────────────────────────────────────────────────────────────────────────────
 # 5. Full-Sample Model Summaries
 #───────────────────────────────────────────────────────────────────────────────
 fit_n_full <- ugarchfit(spec_norm, data=r, solver="hybrid", silent=TRUE)
 fit_s_full <- ugarchfit(spec_t,    data=r, solver="hybrid", silent=TRUE)
 cat("\n--- Full-Sample Normal Model ---\n"); print(fit_n_full)
 cat("\n--- Full-Sample Student's t Model ---\n"); print(fit_s_full)
 
 #───────────────────────────────────────────────────────────────────────────────
 # 6. Plot VaR & ES by Distribution & α
 #───────────────────────────────────────────────────────────────────────────────
 
 # 1) 95% Normal
 p_95_n <- ggplot(df, aes(x = Date)) +
   geom_line(aes(y = Return), color = "grey40", size = 0.3) +
   geom_line(aes(y = VaR5_norm),   color = "red",      size = 0.5) +
   #geom_line(aes(y = ES5_norm),    color = "darkred",  linetype = "dashed", size = 0.8) +
   geom_point(aes(y = vio5_norm), color = "black", shape = 17, size = 2.5, na.rm = TRUE) +
   #geom_point(aes(y = es5_norm),  color = "darkred", shape = 1, size = 2, na.rm = TRUE) +
   labs(title = "VaR at 95% (Normal)", y = "Return / Risk") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5))
 
 # 2) 95% Student’s t
 p_95_t <- ggplot(df, aes(x = Date)) +
   geom_line(aes(y = Return), color = "grey40", size = 0.3) +
   geom_line(aes(y = VaR5_t),    color = "blue",     size = 0.5) +
   #geom_line(aes(y = ES5_t),     color = "darkblue", linetype = "dashed", size = 0.8) +
   geom_point(aes(y = vio5_t), color = "black", shape = 17, size = 2.5, na.rm = TRUE) +
   #geom_point(aes(y = es5_t),  color = "darkblue", shape = 1, size = 2, na.rm = TRUE) +
   labs(title = "VaR at 95% (Student’s t)", y = "Return / Risk") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5))
 
 # 3) 99% Normal
 p_99_n <- ggplot(df, aes(x = Date)) +
   geom_line(aes(y = Return), color = "grey40", size = 0.3) +
   geom_line(aes(y = VaR1_norm),   color = "red",      size = 0.5) +
   #geom_line(aes(y = ES1_norm),    color = "darkred",  linetype = "dashed", size = 0.8) +
   geom_point(aes(y = vio1_norm), color = "black", shape = 17, size = 2.5, na.rm = TRUE) +
   #geom_point(aes(y = es1_norm),  color = "darkred", shape = 1, size = 2, na.rm = TRUE) +
   labs(title = "VaR at 99% (Normal)", y = "Return / Risk") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5))
 
 # 4) 99% Student’s t
 p_99_t <- ggplot(df, aes(x = Date)) +
   geom_line(aes(y = Return), color = "grey40", size = 0.3) +
   geom_line(aes(y = VaR1_t),    color = "blue",     size = 0.5) +
   #geom_line(aes(y = ES1_t),     color = "darkblue", linetype = "dashed", size = 0.8) +
   geom_point(aes(y = vio1_t), color = "black", shape = 17, size = 2.5, na.rm = TRUE) +
   #geom_point(aes(y = es1_t),  color = "darkblue", shape = 1, size = 2, na.rm = TRUE) +
   labs(title = "VaR at 99% (Student’s t)", y = "Return / Risk") +
   theme_minimal() +
   theme(plot.title = element_text(hjust = 0.5))
 
 # Print them one by one:
 print(p_95_n)
 print(p_95_t)
 print(p_99_n)
 print(p_99_t)
 
 # plot function for ES with violations
 plot_es <- function(alpha, dist) {
   # alpha = 0.05 or 0.01
   pct      <- alpha * 100
   es_col   <- paste0("ES",  pct, "_", dist)
   vio_col  <- paste0("es", pct, "_", dist)
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
 # Plot ES
 print(plot_es(0.05, "norm"))  # 95% Normal ES
 print(plot_es(0.05, "t"))     # 95% Student’s t ES
 print(plot_es(0.01, "norm"))  # 99% Normal ES
 print(plot_es(0.01, "t"))     # 99% Student’s t ES
 
 # number of violations
 n_vio_95_n <- sum(!is.na(df$vio5_norm))
 n_vio_95_t <- sum(!is.na(df$vio5_t))
 n_vio_99_n <- sum(!is.na(df$vio1_norm))
 n_vio_99_t <- sum(!is.na(df$vio1_t))
 cat("\nNumber of Violations:\n",
     "95% Normal:", n_vio_95_n, "\n",
     "95% Student’s t:", n_vio_95_t, "\n",
     "99% Normal:", n_vio_99_n, "\n",
     "99% Student’s t:", n_vio_99_t, "\n")
 #number of expected violations
 n_exp_vio_95_n <- nrow(df) * 0.05
 n_exp_vio_95_t <- nrow(df) * 0.05
 n_exp_vio_99_n <- nrow(df) * 0.01
 n_exp_vio_99_t <- nrow(df) * 0.01
 cat("\nExpected Violations:\n",
     "95% Normal:", n_exp_vio_95_n, "\n",
     "95% Student’s t:", n_exp_vio_95_t, "\n",
     "99% Normal:", n_exp_vio_99_n, "\n",
     "99% Student’s t:", n_exp_vio_99_t, "\n")
 # number of ES violations
 n_vio_99_n_es <- sum(!is.na(df$es1_norm))
 n_vio_99_t_es <- sum(!is.na(df$es1_t))
 cat("\nNumber of ES Violations:\n",
     "99% Normal:", n_vio_99_n_es, "\n",
     "99% Student’s t:", n_vio_99_t_es, "\n")
 
 #───────────────────────────────────────────────────────────────────────────────
 # 8. Forecast MSE/RMSE on mean forecasts (Student’s t mean used)
 #───────────────────────────────────────────────────────────────────────────────
 mse <- mean((df$Return - mu_hat)^2, na.rm=TRUE)
 rmse <- sqrt(mse)
 cat("\nMean Forecast MSE:", round(mse,6),
     " RMSE:", round(rmse,6), "\n")
 
 #───────────────────────────────────────────────────────────────────────────────
 # 9. Backtests for both α’s and distributions
 #───────────────────────────────────────────────────────────────────────────────
 vars      <- c("VaR5_norm","VaR5_t","VaR1_norm","VaR1_t")
 alpha_lookup <- c(
   VaR5_norm = 0.05,
   VaR5_t    = 0.05,
   VaR1_norm = 0.01,
   VaR1_t    = 0.01
 )
 res_list  <- vector("list", length(vars))
 
 run_conditional_cc <- function(df_tmp, alpha_target) {
   # df_tmp must have columns I (0/1) and I_lag
   n_hits   <- sum(df_tmp$I)                   # total breaches
   print(n_hits)
   n_trans0 <- sum(df_tmp$I_lag == 0 & df_tmp$I == 1)  # 0->1 transitions
   print(n_trans0)
   n_trans1 <- sum(df_tmp$I_lag == 1 & df_tmp$I == 1)  # 1->1 transitions
   print(n_trans1)
   
   # enough variation?
   if (n_hits == 0 || 
       n_hits == nrow(df_tmp) ||
       n_trans0 == 0 ||
       n_trans1 == 0) {
     return(c(chisq = NA, p = NA))
   }
   
   # fit the regression
   mod_ind <- lm(I ~ 1 + I_lag, data = df_tmp)
   
   # compute robust covariance
   V <- vcovHC(mod_ind, "HC1")
   b <- coef(mod_ind)
   
   # joint Wald test H0: intercept = alpha_target; slope = 0
   w <- wald.test(b = b, Sigma = V,
                  Terms = 1:2,  # indices of (Intercept) and I_lag
                  H0    = c(alpha_target, 0))
   
   return(c(chisq = w$result$chi2["chi2"],
            p     = w$result$chi2["P"]))
 }
 # Helper for the LR‐based conditional‐coverage test
 lr_conditional_cc <- function(hits, alpha) {
   hits1 <- hits[-length(hits)]
   hits2 <- hits[-1]
   n00 <- sum(hits1==0 & hits2==0)
   n01 <- sum(hits1==0 & hits2==1)
   n10 <- sum(hits1==1 & hits2==0)
   n11 <- sum(hits1==1 & hits2==1)
   # Unrestricted MLEs
   pi0_hat <- if((n00+n01)>0) n01/(n00+n01) else NA
   pi1_hat <- if((n10+n11)>0) n11/(n10+n11) else NA
   ll_unres <- n00*log(1-pi0_hat) + n01*log(pi0_hat) +
     n10*log(1-pi1_hat) + n11*log(pi1_hat)
   # Restricted at alpha
   ll_res   <- (n00+n10)*log(1-alpha) + (n01+n11)*log(alpha)
   LR        <- 2*(ll_unres - ll_res)
   p         <- 1 - pchisq(LR, df=2)
   c(LR, p)
 }
 
 for(i in seq_along(vars)) {
   var_col      <- vars[i]
   alpha_target <- alpha_lookup[[var_col]]
   
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
   
   cc_res <- run_conditional_cc(df_tmp, alpha_target = alpha_target)
   print(var_col)
   print(cc_res)
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
     # use HC1 robust covariance
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
 #print(results_df)
 
 
 # ──────────────────────────────────────────────────────────────────────────────
 # 8. Residual Diagnostics
 # ──────────────────────────────────────────────────────────────────────────────
 
 # 1. Extract standardized residuals
 res_n <- residuals(fit_n_full, standardize = TRUE)
 res_t <- residuals(fit_s_full,    standardize = TRUE)
 nu    <- coef(fit_s_full)["shape"]  # estimated degrees of freedom
 
 # 2. Histogram + density for Normal model
 df_n <- data.frame(res = res_n)
 p_hist_n <- ggplot(df_n, aes(x = res)) +
   geom_histogram(aes(y = ..density..), bins = 50, fill = "grey80", color = "black") +
   stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                 color = "red", size = 1) +
   labs(title = "Histogram of Std. Residuals (Normal GARCH)",
        x = "Standardized Residual", y = "Density") +
   theme_minimal()
 
 # 3. Histogram + density for t model
 df_t <- data.frame(res = res_t)
 p_hist_t <- ggplot(df_t, aes(x = res)) +
   geom_histogram(aes(y = ..density..), bins = 50, fill = "grey80", color = "black") +
   stat_function(fun = dt, args = list(df = nu),
                 color = "blue", size = 1) +
   labs(title = paste0("Histogram of Std. Residuals (std-GARCH, ν=", round(nu,2), ")"),
        x = "Standardized Residual", y = "Density") +
   theme_minimal()
 
 # 4. QQ-plot for Normal residuals
 p_qq_n <- ggplot(df_n, aes(sample = res)) +
   stat_qq(color = "black") +
   stat_qq_line(color = "red") +
   labs(title = "QQ Plot: Normal GARCH Residuals",
        x = "Theoretical Quantiles", y = "Sample Quantiles") +
   theme_minimal()
 
 # 5. QQ-plot for t residuals
 p_qq_t <- ggplot(df_t, aes(sample = res)) +
   stat_qq(distribution = qt, dparams = list(df = nu), color = "black") +
   stat_qq_line(distribution = qt, dparams = list(df = nu), color = "blue") +
   labs(title = paste0("QQ Plot: std-GARCH Residuals (ν=", round(nu,2), ")"),
        x = "Theoretical Quantiles", y = "Sample Quantiles") +
   theme_minimal()
 
 # 6. Arrange all four plots in a 2×2 grid
 grid.arrange(p_hist_n, p_hist_t, p_qq_n, p_qq_t, ncol = 2, nrow = 2)
 