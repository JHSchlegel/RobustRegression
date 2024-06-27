################################################################################
###                       Robust Regression                                  ###
################################################################################





























################################################################################
###                       Non-Linear Regression                              ###
################################################################################


# Confidence Intervals ----------------------------------------------------
## Linearization:
sm.nls = summary(fit.nls)
h = qt(0.975, df = nrow(chlor) - 3) * sm.nls$coefficients[, 2] 
ci = coef(fit.nls) + cbind(-h, h)

## Profiling:
confint(fit.nls)



# Residual Analysis -------------------------------------------------------
## TA-plot:
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


## QQ-Plot:
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



h = qqnorm(resid(fit.nls))
qqline(resid(fit.nls))
identify(h)



