################################################################################
###                       Robust Regression                                  ###
################################################################################


# M estimation ------------------------------------------------------------
## HuberM:
robfit = huberM(df, se = T, k = 1.345)
robfit$mu + c(-1, 1) * qt(0.975, length(df) - 1) * robfit$SE

## rlm:
rlm(ValuesOrg ~ Block + Variety, data = oats, psi = psi.huber, method = "M", maxit = 50)

# MM estimation -----------------------------------------------------------
fit.mm = lmrob(formula, data = dat)


# SMDM Estimation ---------------------------------------------------------
lmrob(y~., data = synth, setting = "KS2014")



# ANOVA -------------------------------------------------------------------
## Robust LM:
anova(m1, m2, test = "Wald")
anova(m1, m2, test = "Deviance")

## Robust GLM:
anova(ep.rglm, ep.rglm2, test = "QD")


# Variable Selection ------------------------------------------------------
## lm:
step(fit.lm, direction= "backward", trace = F)

## roblm:
h.cont = lmrobdet.control(bb=0.5, efficiency = 0.85, family = "bisquare")
rlm1 = lmrobdetMM(SpGew~., data = wood.cont, control = h.cont)
step.lmrobdetMM(rlm1)




# Robust GLM --------------------------------------------------------------
## No potential outliers in x-space:
s.rglm = glmrob(PASS ~., data = apt, family = binomial(), method = "Mqle")
plot(s.rglm$w.r); plot(s.rglm$w.x)

## Potential outliers in x-space:
s.rglm = glmrob(ep.formula, data = epilepsy, family = poisson(), method = "Mqle", weights.on.x = "hat")
plot(s.rglm$w.r); plot(s.rglm$w.x)


## Plotting w/ obs. number:
ep.df = data.frame(wx = ep.rglm2$w.x, label = 1:length(ep.rglm2$w.x))

ep.df |> 
  ggplot(aes(x=label, y = wx)) +
  geom_point() +
  geom_text(aes(label = label), hjust = 1.5, vjust = 0.5) +
  theme_bw()



# Robust Covariance -------------------------------------------------------
(cov.rob <- CovRobust(d.wood, control = "mcd"))


# Mahalanobis Distance ----------------------------------------------------
## classical:
cov.class = cov(d.wood)
md.class = mahalanobis(d.wood, center = colMeans(d.wood), cov = cov.class)

##  robust:
cov.rob = CovRobust(d.wood, control = "mcd")
md.rob = mahalanobis(d.wood, center = cov.rob@center, cov = cov.rob@cov)

## Plotting:
# robust vs classical:
df.md |> 
  pivot_longer(-obs, names_to = "type") |> 
  ggplot(aes(x = obs, y = value, col=type, shape = type)) +
  geom_point() +
  theme_bw()

# Chisq. QQ-Plot:
df.md |> 
  mutate(label = 1:nrow(df.md)) |> 
  arrange(robust) |> 
  mutate(theoretical_quantiles = qchisq(ppoints(nrow(df.md)), df = 5)) |> 
  ggplot(aes(x = theoretical_quantiles, y = robust)) +
  geom_point() + 
  geom_text(aes(label = label), vjust = -0.5, hjust = 0.5) +
  ylim(0, 250) +
  geom_abline(aes(intercept = 0,slope =1), linewidth = 0.2, color = "red") +
  theme_bw() +
  labs(
    title = "Robust MCD Covariance Estimate",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  )



# PCA ---------------------------------------------------------------------
## Classical:
pca.class = princomp(d.wood, cor = T)

## Robustly standardized:
stand_rob_loc_scale = scale(x = d.wood, center = colMedians(as.matrix(d.wood)), scale = apply(d.wood, 2, mad))
pca.rob_loc_scale = princomp(stand_rob_loc_scale)

## Robustly estimated covmat:
pca.mcd = rrcov::PcaCov(d.wood, scale = T)

## get loadings:
loadings(pca.obj)

## pairs plots:
pairs(pca.obj$scores)

## Scree plots:
# base R:
screeplot(pca.obj)

# ggplot:
pca.var = data.frame(
  component = 1:5,
  classical = pca.class$sdev^2,
  robust_scaling = pca.rob_loc_scale$sdev^2,
  mcd = pca.mcd$eigenvalues
)

pca.var |> 
  pivot_longer(-component, names_to="type") |> 
  ggplot(aes(x = component, y = value, color = type, linetype = type, shape = type)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(
    x = "Component",
    y = "Variance"
  )



# LDA ---------------------------------------------------------------------
## Classical:
lda.class = MASS::lda(Klasse ~ ., data = klasse.df)

## Only W robustly:
lda.mve = MASS::lda(Klasse~ ., data = klasse.df, method = "mve")

## W and locations of groups robustly (see utils file):
rlda.groups = rlda(x = klasse.df[, 2:3], grouping = klasse.df$Klasse)

## Scatter plot first and secondend discriminant axis
p.ldv(rlda.groups, data = klasse.df[, 2:3], group = klasse.df$Klasse)
title("RLDA")

## Plot decision boundary
p.predplot(rlda.groups, data = klasse.df[, 2:3], group = klasse.df$Klasse)
title("RLDA")

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
## Scatter plots by variable:
scatter.smooth(dat$var1, resid(fit))

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
# percentile:
summary(boot.nls); boot.nls$bootCI

# Wald:
confint2(boot.nls, method = "asymptotic")

## Profile Likelihood:
confint(boot.nls)


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



# Boot Library ------------------------------------------------------------

library(boot)
f.bod <- function(rs, ind){
  bsY <- fitted(D.bod.nls) + rs[ind]
  coef(nls(bsY ~ Th1*(1-exp(-Th2*days)), data=D.bod,
           start=coef(D.bod.nls)))
}
h.rs <- scale(resid(D.bod.nls), scale=FALSE) ## mean centred residuals
## bootstrapping
set.seed(seed=117)
D.bod.nls.Boot2 <- boot(h.rs, f.bod, R=4999, stype="i" )
## bias-corrected accelerated (bca) bootstrap interval
rbind(boot.ci(D.bod.nls.Boot2, conf=0.95, type="bca", index=1)$bca[4:5],
       boot.ci(D.bod.nls.Boot2, conf=0.95, type="bca", index=2)$bca[4:5])





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
plot(calibr.int)

## nonparametric bootrstrap:
calibr.int = invest(
  fit.nls, y0 = 5.1, interval = "percentile", 
  boot.type ="nonparametric", nsim = 300, seed = 42
)
plot(calibr.int)



