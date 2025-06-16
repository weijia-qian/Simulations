library(splines)
library(tidyverse)

s <- seq(0, 1, length.out = 500)
B <- bs(s, df = 8, intercept = TRUE)
bs_coef = c(0.2, 0.4, 0.2, 0.2, 0.2, 0.6, 0.3, 0.2)
beta_1 <- B %*% bs_coef
tibble(s = s,
       beta_1 = beta_1) %>%
  ggplot(aes(s, beta_1)) +
  geom_line()

load('/Users/weijia/Research/FDA/Simulations/Output/20241227/10.RDA')
df_coef <- bind_rows(lapply(results, function(x) x$coef)) %>%

  df_coef <- data.frame()

  
  for (iteration in results){

    
    # Extract data for the current iteration
    time <- seq(0, 1, length.out = iteration$info$nS)
    tmp_data <- data.frame(scenario = iteration$info$scenario,
                           iter = iteration$info$iter,
                           family = iteration$info$family,
                           n = iteration$info$n,
                           nS = iteration$info$nS,
                           censor_rate = iteration$info$censor_rate,
                           time = time,
                           beta_type = iteration$info$beta_type,
                           true_coef = iteration$coef$true_coef,
                           est_coef_norm = iteration$coef$est_coef_norm,
                           # lb_coef_norm = iteration$coef$lb_coef_norm,
                           # ub_coef_norm = iteration$coef$ub_coef_norm,
                           # cma_lb_coef_norm = iteration$coef$cma_lb_coef_norm,
                           # cma_ub_coef_norm = iteration$coef$cma_ub_coef_norm,
                           est_coef_cox = iteration$coef$est_coef_cox)
                           # lb_coef_cox = iteration$coef$lb_coef_cox,
                           # ub_coef_cox = iteration$coef$ub_coef_cox,
                           # cma_lb_coef_cox = iteration$coef$cma_lb_coef_cox,
                           # cma_ub_coef_cox = iteration$coef$cma_ub_coef_cox)
    
    # Append to all_coef
    df_coef <- rbind(df_coef, tmp_data)

  }


save(df_coef_mise, df_coef, file = "coef.RDA")
