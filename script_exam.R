################################################################################
###                       Robust Regression                                  ###
################################################################################





























################################################################################
###                       Non-Linear Regression                              ###
################################################################################


# Summary -----------------------------------------------------------------
summary(fit.nls, corr = T)


# Confidence Intervals ----------------------------------------------------
## Linearization:
sm.nls = summary(fit.nls)
h = qt(0.975, df = nrow(chlor) - 3) * sm.nls$coefficients[, 2] 
ci = coef(fit.nls) + cbind(-h, h)

## Profiling:
confint(fit.nls)



# Residual Analysis -------------------------------------------------------
## Tukey-Anscombe plot:
# should not show any special structure

# ggplot:
data.frame(
  fitted = fitted(fit.nls), 
  resid = resid(fit.nls),
  label = 1:nrow(chlor)
) |> 
  ggplot(aes(x = fitted, y = resid)) + 
  geom_point(size = 0.5) +
  geom_smooth(color = "firebrick", fill = "firebrick", alpha = 0.5, linewidth =0.5) +
  geom_hline(aes(yintercept=0), color = "grey40", linewidth = 0.5, linetype = "dashed")+
  geom_text(aes(label = label), hjust = 1.5, vjust = 0.5, size = 2) +
  theme_bw()

# base R:
plot(fitted(fit.nls), resid(fit.nls))
abline(h=0, lty=2)

## QQ-Plot:
#ggplot:
data = data.frame(resid = resid(fit.nls), label = 1:nrow(chlor)) |> 
  arrange(resid) |>
  mutate(
    theoretical_quantiles = qnorm(ppoints(nrow(chlor))),
    sample_quantiles = resid
  ) 

slope <- sd(data$sample_quantiles) / sd(data$theoretical_quantiles)
intercept <- mean(data$sample_quantiles) - slope * mean(data$theoretical_quantiles)

data |> 
  ggplot(aes(x = theoretical_quantiles, y = sample_quantiles)) +
  geom_point(size = 0.5) +
  geom_text(aes(label = label), size = 2, hjust = 1.5, vjust = 0.5) +
  geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") +
  theme_bw()


# base R:
h = qqnorm(resid(fit.nls))
qqline(resid(fit.nls))
identify(h)



# Likelihood Profile Traces -----------------------------------------------
r.prof = profile(fit.nls)
p.profileTraces(r.prof)
# diagonal elements: profile t-plots
# likelihood profile traces: lower diagonal



# Bootstrapping -----------------------------------------------------------
boot.nls = nlsBoot(fit.nls, niter = 999)

## CI:
summary(boot.nls); boot.nls$bootCI


## Manually (w/o quantiles):
# Sort the bootstrap samples
sorted_bootstrap_samples <- sort(bootstrap_samples)

# Calculate the indices for the 2.5th and 97.5th percentiles
lower_index <- ceiling(0.025 * n_boot)
upper_index <- floor(0.975 * n_boot)

# Extract the confidence interval values
ci_lower_sorted <- sorted_bootstrap_samples[lower_index]
ci_upper_sorted <- sorted_bootstrap_samples[upper_index]


## Pairs plot:
plot(boot.nls)


## Histogram:
data.frame(boot.nls$coefboot, iter = 1:991) |> 
  pivot_longer(-iter, names_to = "param") |> 
  ggplot(aes(x = value)) +
  geom_histogram(fill = "midnightblue", color = "white") +
  facet_wrap(~param, scale = "free_x") +
  theme_bw()



# Prediction Intervals ----------------------------------------------------
newdat = list(dose = c(-1.5, 0, 2))
conf.int = investr::predFit(fit.nls, newdata = newdat, interval = "confidence")
preds.int = investr::predFit(fit.nls, newdata = newdat, interval = "prediction")
# width:
preds.int[, 3] - preds.int[2]



# Calibration Intervals ---------------------------------------------------
## Wald approach:
# measured a single value:
calibr.int = invest(fit.nls, y0 = 5.1, interval = "inversion")
# multiple measurements fo rthe same lot -> smaller int:
calibr.int = invest(fit.nls, y0 = c(5.1, 7.3, 6.2), interval="inversion")


## bca bootstrap:
calibr.int = invest(
  fit.nls, y0 = 5.1, interval = "percentile", 
  boot.type ="parametric", nsim = 300, seed = 42
)

## nonparametric bootrstrap:
calibr.int = invest(
  fit.nls, y0 = 5.1, interval = "percentile", 
  boot.type ="nonparametric", nsim = 300, seed = 42
)




