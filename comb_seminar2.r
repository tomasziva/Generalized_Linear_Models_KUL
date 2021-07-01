rm(list = ls())

# packages needed
# install.packages("countreg", repos="http://R-Forge.R-project.org")
pkg <- c("ggplot2", "dplyr", "AER", "DHARMa", "sandwich", "lmtest", "countreg", "pscl")
invisible(sapply(pkg, library, character.only = TRUE))

# read the data
victim <- read.table("victim.txt", header = TRUE)
#str(victim)

# frequency table
freq <- cbind(table(victim$race, victim$resp), "ttl" = table(victim$race))
freq

# relevel the race 
victim$race <- relevel(victim$race, ref = "white")

# Poisson model -----------------------------------------------
mod <- glm(resp~race, family = poisson, data = victim)
summary(mod)

# risk ratio and CI
exp(cbind(RR = mod$coefficients, confint(mod)))

# the predicted mean responses
pred <- predict(mod, newdata = data.frame(race = c("white", "black")),
                type = "response", se.fit = TRUE)

# fitted counts using the model mean estimates
fitted_w <- dpois(0:6, lambda = pred$fit[[1]])*freq["white", "ttl"]
fitted_b <- dpois(0:6, lambda = pred$fit[[2]])*freq["black", "ttl"]
freq_poisson <- cbind(obs_w = freq["white", 1:7],  fitted_w = round(fitted_w, 2), 
                      obs_b = freq["black", 1:7], fitted_b = round(fitted_b, 2))
freq_poisson

# rootogram
countreg::rootogram(mod)

# significance test
# likelihood ratio test
anova(mod, test = "LRT")
# Wald test
car::Anova(mod, test = "Wald", type = 3)
# score test
anova(mod, test = "Rao")

# GOF
# Pearson test
x2 <- sum(residuals(mod, type = "pearson")^2)
n <- nrow(victim)
p <- length(coef(mod))
data.frame(X2 = x2, df = n - p, pvalue = (1-pchisq(x2, n-p)))

# Deviance test
data.frame(deviance = mod$deviance, df = mod$df.residual, pvalue = (1-pchisq(mod$deviance, mod$df.residual)))

# Residual analysis
r_dev <- residuals(mod, type = "deviance")
plot(victim$race, r_dev, xlab = "Race", ylab = "Deviance Residuals")
abline(h=0,lty=2,col="grey")


# overdispersion ----------------------------------------------
# raw variance
list(obs_mean = tapply(victim$resp, victim$race, mean), obs_variance = tapply(victim$resp, victim$race, var))

# dispersion test
dispersiontest(mod)
sim <- simulateResiduals(mod, refit = TRUE)
testDispersion(sim)

# sandwich estimator
coeftest(mod, vcov. = sandwich, test ="Chisq")

# quasi-likelihood
mod_quasi <- glm(resp~race, family = quasipoisson, data = victim)
summary(mod_quasi)

# negative binomial
mod_nb <- glm.nb(resp~race, data = victim)
summary(mod_nb)

# Pearson
x2 <- sum(residuals(mod_nb, type = "pearson")^2)
df = mod_nb$df.residual
data.frame(X2 = x2, df = df, pvalue = pchisq(x2, df, lower.tail = FALSE))
# Deviance test
dev <- summary(mod_nb)$deviance
data.frame(dev = dev, df = df, pvalue = pchisq(dev, df, lower.tail = FALSE))

# likelihood ratio (nb vs poisson)
lrtest(mod, mod_nb)

# rootogram
countreg::rootogram(mod_nb)

# quasi-likelihood vs negative binomial
mean <- pred$fit
disp <- cbind(race = c("white", "black"),
              round(data.frame(mean = mean, var_obs = with(victim, tapply(resp, race, var)),
                               var_ql = summary(mod_quasi)$dispersion * mean,
                               var_nb = mean + (1/mod_nb$theta)*mean^2), 3))
# plot mean-variance 
with(disp, plot(mean, var_obs, ylim = c(0, 2), pch = 19))
with(disp, points(mean, var_ql, pch = 8))
with(disp, points(mean, var_nb, pch = 7))
legend("topleft", pch = c(19, 8, 7), legend = c("Observed", "Quasi-likelihood", "Negative binomial"))

# summarize all the results
round(data.frame(po = coef(mod),  po_sw = coef(mod),  po_ql = coef(mod_quasi),
                 po_nb = coef(mod_nb),  se.po = summary(mod)$coefficients[, 2],
                 se.po_sw = coeftest(mod, vcov. = sandwich)[, 2],
                 se.po_ql = summary(mod_quasi)$coefficients[, 2],
                 se.po_nb = summary(mod_nb)$coefficients[, 2] ), 3)

# zero inflation -------------------------------------------------------
testZeroInflation(sim)
freq_poisson[1,]

# prob(0) ~ 1
mod_zip <- zeroinfl(resp ~ race | 1, data = victim)
summary(mod_zip)

# prob(0) ~ race
mod_zip2 <- zeroinfl(resp ~ race | race, data = victim)
# negative binomial + prob(0) ~ race
mod_zinb2 <- zeroinfl(resp ~ race | race, data = victim, dist = "negbin")
summary(mod_zip2)

cbind(AIC(mod, mod_nb, mod_zip, mod_zip2, mod_zinb2), 
      BIC(mod, mod_nb, mod_zip, mod_zip2, mod_zinb2))
