# packages needed
pkg <- c("ggplot2", "dplyr", "gridExtra", "car", "multcomp", "dunn.test", "gplots", "brglm2",
         "caret", "BaylorEdPsych", "AICcmodavg", "pROC", "tidyr")
invisible(sapply(pkg, library, character.only = TRUE))

# Read data
income <- read.table("income.txt", header = TRUE, sep = "", dec = ".")

# standardize age (agest) and hours (hourst) 
incomest <- income %>% dplyr::select(-age, -hours)
incomest$agest <- as.numeric(scale(income$age))
incomest$hourst <- as.numeric(scale(income$hours))

# Plot gender versus income
incomest %>%
  group_by(income, gender) %>%
  tally() %>%
  group_by(income) %>%
  mutate(proportion = n / sum(n)) %>%
  ggplot() + geom_col(aes(x = factor(gender),  y = proportion, fill = factor(income)), position = "fill")

ggplot(incomest,  aes(x = gender,  fill = income)) +  geom_bar(position = "stack")


# Plot race versus income
incomest %>%
  group_by(income, race) %>%
  tally() %>%
  group_by(income) %>%
  mutate(proportion = n / sum(n)) %>%
  ggplot() + geom_col(aes(x = factor(race),  y = proportion, fill = factor(income)), position = "fill")

ggplot(incomest, aes(x = race, fill = income)) +  geom_bar(position = "stack")

# gender versus working time
ggplot(incomest, aes(y=hourst, fill = income)) + geom_boxplot()
ggplot(incomest, aes(hourst)) + geom_histogram() + facet_grid(gender~.)

## Plot distribution working time by education
ggplot(incomest, aes(x=hourst)) +  geom_density(alpha=.2, fill="purple4") + facet_wrap( ~ edu, ncol = 3)


ggplot(incomest, aes(edu, hourst, fill = edu)) + geom_boxplot() + theme(legend.position = "none")

# work class vs income
ggplot(incomest, aes(workclass, fill = income)) + geom_bar(position = "fill")
# marital status vs income
ggplot(incomest, aes(mari.sta, fill = income)) + geom_bar(position = "fill")

# ANOVA working hours vs edu -------------------------------------------------
incomest$edu<-factor(incomest$edu, levels = c("dropout", "HighGrad", "Bachelors", "Community", "Master", "PhD"))
model.fit<-lm(hourst~edu, data=incomest)
summary(model.fit)
pair <- glht(model.fit, linfct = mcp(edu = "Tukey"))
par(mar = c(4, 10, 4, 2)+0.1)
plot(confint(pair))

# Normality check
ks.test(residuals(model.fit), "pnorm", alternative = "two.sided")
qqnorm(residuals(model.fit))
qqline(residuals(model.fit))
# homoscedasticity check
leveneTest(hourst~edu, data = incomest)
# non-parametric test Kruscal-Wallis test
kruskal.test(hourst ~ edu, data = incomest)
dunn.test(incomest$hourst, incomest$edu, method = "bonferroni")

## Create a train data set containing 80% of the data. ----------------------
set.seed(8)
train_index <- sample(1:nrow(incomest), 0.8 * nrow(incomest))
data.train <- incomest[train_index,]
## Create a test data set containing 20% of the data.
test_index <- setdiff(1:nrow(incomest), train_index)
data.test <- incomest[test_index,]

# logit model ----------------------------------------
data.train$race<-relevel(data.train$race,ref="White")
data.train$edu<-relevel(data.train$edu,ref="dropout")

# quasi-complete separation detection
glm(income ~ race + edu + mari.sta + workclass + gender + hourst + agest, 
    data=data.train,family=binomial(link="logit"), method = "detect_separation",
    linear_program = "dual")

# linear model
# PENALIZED ML
m.pred<-glm(income ~ race + edu + mari.sta + workclass + gender + hourst + agest,
            data=data.train,family=binomial(link="logit"),
            method = "brglmFit", maxit = 1000)
## NOTE! The model with penalized ML cannot converge. 
## So we had to switch back to ML for the rest of the project.
m.pred<-glm(income ~ race + edu + mari.sta + workclass + gender + hourst + agest,
            data=data.train,family=binomial(link="logit"))
summary(m.pred)
Anova(m.pred, test = "LR", type = 3)
Anova(m.pred, test = "Wald", type = 3)

# GOF
# Pearson
x2 <- sum(residuals(m.pred, type = "pearson")^2)
df <- m.pred$df.residual
data.frame(X2 = x2, df = df, pvalue = 1 - pchisq(x2, df))
# deviance
dev <- m.pred$deviance
data.frame(deviance = dev, df = df, pvalue = 1 - pchisq(dev, df))


# quadratic model 
m.pred2<-glm(income ~ race + edu + mari.sta + workclass + gender + poly(hourst, 2) + poly(agest, 2),  data=data.train,family=binomial(link="logit"))
summary(m.pred2)
# LRT
anova(m.pred, m.pred2, test = "Chisq")
anova(m.pred, m.pred2,test='Rao')


# prediction performance ----------------------------------------
# linear model
data.test$greP<-predict.glm(m.pred,newdata=data.test,type="response")
data.test$pred <- ifelse(data.test$greP >=0.25, ">50K", "<=50K")
# classification table
table(data.test$income, data.test$pred)

# quadratic model
data.test$greP2<-predict(m.pred2,newdata=data.test,type="response")
data.test$pred2 <- ifelse(data.test$greP2 >=0.5, ">50K", "<=50K")

conf_matrix <- table(data.test$pred, data.test$income)
specificity(conf_matrix) #true negative rate
sensitivity(conf_matrix) #true positive rate
# accuracy
dt.tbl <- as.data.frame(table(data.test$income, data.test$pred))
sum(dt.tbl$Freq[dt.tbl$Var1 == dt.tbl$Var2]) / sum(dt.tbl$Freq)

# quadratic model 
data.test$greP2<-predict(m.pred2,newdata=data.test,type="response")
data.test$pred2 <- ifelse(data.test$greP2 >=0.5, ">50K", "<=50K")

dt.melt2 <- data.test %>% select(income, pred2) %>% gather()
table(dt.melt2$value, dt.melt2$key)
conf_matrix2 <- table(data.test$pred2, data.test$income)
specificity(conf_matrix2) #true negative rate
sensitivity(conf_matrix2) #true positive rate

# roc
roccurve <- roc(data.test$income, data.test$greP)
roccurve2 <- roc(data.test$income, data.test$greP2)
plot(roccurve, col = "green")
plot(roccurve2, add = TRUE, col = "red")

# optimal cut-off (with the linear model)
cutoff <- cutpointr(data.test, pred, income,  direction = ">=",   pos_class = 1,  neg_class = 0, 
                    method = maximize_metric,  metric = sum_sens_spec)
summary(cutoff)
plot(cutoff)

# prediction comparison
# aic
aic <- AIC(m.pred, m.pred2)
# pseudo R squred
McFadden <- c(BaylorEdPsych::PseudoR2(m.pred)[1],
              BaylorEdPsych::PseudoR2(m.pred2)[1])
Cox.Snell <- c(BaylorEdPsych::PseudoR2(m.pred)[3],
               BaylorEdPsych::PseudoR2(m.pred2)[3])
Nagelkerke <- c(BaylorEdPsych::PseudoR2(m.pred)[4],
                BaylorEdPsych::PseudoR2(m.pred2)[4])
cbind(aic, McFadden, Cox.Snell, Nagelkerke)

# Akaike weights
m.pred.4 <- glm(income ~ agest + hourst + edu + gender, family = binomial(link = "logit"),
                data = data.train)
m.pred.x <- glm(income ~ agest + hourst + edu + race + mari.sta + workclass + gender +
                  edu:gender + mari.sta:gender, family = binomial(link = "logit"),
                data = data.train)
mod_list <- list(m.pred.4, m.pred, m.pred.x, m.pred2)
mod_names <- c("Four", "Linear", "Interaction", "Quadratic")

# Akaike weights
AICcmodavg::aictab(cand.set = mod_list, modnames = mod_names)
