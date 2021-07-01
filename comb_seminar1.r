library(aod)
library(car)
library(ggplot2)
library(dplyr)

djasem <- read.table("djaSem.txt", header = TRUE)

# calculate proportions
djasem$prop <- djasem$y / djasem$n  
boxplot(prop~group, data = djasem)
points(1:2, tapply(djasem$prop, djasem$group, mean), pch = 16, col = "red")

colnames(djasem)[1] <- "treatment"

# logit model --------------------------------------------------------
mod_bn <- glm(cbind(y, n-y)~treatment, family = binomial(link = "logit"), data = djasem)
summary(mod_bn)

# odds ratios
exp(mod_bn$coefficients)
# combined 
exp(cbind(OR = mod_bn$coefficients, confint(mod_bn)))

# likelihood test
anova(mod_bn, test = "Chisq")
# Wald test
Anova(mod_bn, test = "Wald", type = 3)
# score test
anova(mod_bn, test = "Rao")

# predicted prob. death
predict(mod_bn, newdata = data.frame(treatment = c("CTRL", "TREAT")), se.fit = TRUE,
        type = "response")

# plot pred. prob. of death
pred <- predict(mod_bn, newdata = data.frame(treatment = c("CTRL", "TREAT")), se.fit = TRUE)  
pred_df <- data.frame(treatment = c("CTRL", "TREAT"), fit = pred$fit, se = pred$se.fit)

ggplot(pred_df, aes(treatment, exp(fit)/(1+exp(fit)))) + geom_point() + 
  geom_errorbar(aes(ymin = exp(fit-1.96*se)/(1+exp(fit-1.96*se)), 
                    ymax = exp(fit+1.96*se)/(1+exp(fit+1.96*se)),width = 0.05)) +
  labs(x = "Group", y = "Predicted Probability of Death (95% CI)")

# GOF
# Pearson
x2 <- sum(residuals(mod_bn, type = "pearson")^2)
data.frame(X2 = x2, df = mod_bn$df.residual, pvalue = 1 - pchisq(x2, 75-2))

# deviance
data.frame(deviance = mod_bn$deviance, df = mod_bn$df.residual,  pvalue = 1-pchisq(mod_bn$deviance, mod_bn$df.residual))

model.rdev = residuals(mod_bn, type = "deviance")
plot(predict(mod_bn), model.rdev, xlab = "Linear predictor", ylab = "Deviance Residuals")
abline(h=0,lty=2,col="grey")


loess.model.rdev <- loess(model.rdev ~ predict(mod_bn))
model.lo.pred.dev <- predict(loess.model.rdev, se=T)
j.model = order(predict(mod_bn))
lines(predict(mod_bn)[j.model], model.lo.pred.dev$fit[j.model], col = "blue", lwd = 3)
lines(predict(mod_bn)[j.model], model.lo.pred.dev$fit[j.model] 
      +2*model.lo.pred.dev$s[j.model], lty = 2, col = "red")
lines(predict(mod_bn)[j.model], model.lo.pred.dev$fit[j.model]
      -2*model.lo.pred.dev$s[j.model], lty = 2, col = "red")


# Delta Deviance, Chi-Sq & Cook's Distance Plots~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model.hat = hatvalues(mod_bn)
model.rs = rstudent(mod_bn)
model.st = residuals(mod_bn, type="pearson")/sqrt(1 - model.hat)
model.rs.pearson = residuals(mod_bn, type = "pearson")

deltaX2 = model.st^2
deltaD = model.hat*model.st^2+residuals(mod_bn, type="deviance")^2
cook.dist = (model.rs.pearson^2)*model.hat/(2*(1-model.hat)^2)

# Delta Chi-squared    Outliers: 1 34 39 71  (Change title if you change the model)
plot(deltaX2, xlab="Case number", ylab="Delta Chi-Square", main = 'Delta Chi-Square (Logit)', col = 'red', pch=20)
lines(deltaX2)
identify(deltaX2)

# Delta Deviance      Outliers: 1 34 71
plot(deltaD, xlab="Case number", ylab="Delta Deviance", main = 'Delta Deviance (Logit)',col = 'red', pch=20)
lines(deltaD)
# identify(deltaD)

# Cook Distance      Outliers: 1 34 39 71
plot(cook.dist, xlab="Case number", ylab="Cook's Distance", main = "Cook's Distance (Logit)", col = 'red', pch=20)
lines(cook.dist)
# identify(cook.dist)

# Look at influential Points~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# influential observations
influ.obs = djasem[c(1, 34, 39, 71), ]

# Avg # deaths (grp vs influential obs)
influ.obs$y   #influential obs
# aggregate(djasem$y, by=list(Category=djasem$group), FUN=mean)  # mean num deaths by grp
djasem %>% group_by(treatment) %>% summarise(count = mean(y))

# Avg # no deaths (grp vs influential obs)
influ.obs$n - influ.obs$y   #influential obs
djasem$ny = djasem$n - djasem$y
# aggregate(djasem$ny, by=list(Category=djasem$group), FUN=mean)  # mean num no deaths by grp
djasem %>% group_by(treatment) %>% summarise(count = mean(n))


# quasi-likelihood ---------------------------------------
mod_ql <- glm(cbind(y, n-y)~treatment, family = quasibinomial(link = "logit"), data = djasem)
summary(mod_ql)

# Pearson
x2 <- sum(residuals(mod_ql, type = "pearson")^2)/1.294293
df <- summary(mod_ql)$df.residual
data.frame(X2 = x2, df = df, pvalue = 1 - pchisq(x2, df))
# deviance
deviance <- mod_ql$deviance/1.294293
data.frame(deviance = deviance, df = mod_ql$df.residual, pvalue = 1 - pchisq(deviance, df))
# LRT
anova(mod_bn, mod_ql, test = "LRT")

# beta-binomial
mod_beta <- betabin(cbind(y, n - y) ~ treatment, random = ~ 1, data = djasem)
mod_beta

# predicted probability
predict(mod_beta, newdata = data.frame(treatment = c("CTRL", "TREAT")), type = "response")

# compare AIC
# AIC(mod_bn, mod_beta)
