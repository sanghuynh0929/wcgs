# Refs: https://cpb-us-w2.wpmucdn.com/voices.uchicago.edu/dist/6/1366/files/2017/01/Statistics-Day4.pdf
# https://thomaselove.github.io/431-notes/17-wcgs.html
# https://thomaselove.github.io/431-notes/31-prediction_wcgs.html 


# Load required libraries


# install.packages("epitools")
library(epitools)

# Load the WCGS dataset
data(wcgs)

# install.packages("psych")
library(psych)
library(ggplot2)
library(tidyr)
library(faraway)
library(dplyr)
library(car)
library(survival)
library(survminer)
library(Hmisc)
library(bestglm)
library(lmtest)

# View the dataset
View(wcgs)
attach(wcgs)


wcgs$smoker <- ifelse(wcgs$ncigs0 > 0, 1, 0)
wcgs$bmi <- 703 * wcgs$weight0 / (wcgs$height0) ^ 2

pairs.panels(wcgs[, -1])

# chol co outlier
# sbp, dbp co correlation



# Suppose that after consulting with clinical experts, we want to ensure that:

# all age values are between 39 and 59 years
# all bmi are between 15 and 50
# all sbp are between 80 and 250 mm Hg
# all dbp are between 50 and 200 mm Hg
# all values of sbp-dbp are at least 10 and no more than 90 mm Hg
# all values of chol are between 100 and 400 mg/dl


Hmisc::describe(wcgs %>% select(age0, bmi, sbp0, dbp0, ncigs0, smoker, dibpat0, chol0))
# wcgs %>% select(age0, bmi, sbp0, dbp0, ncigs0, smoker, dibpat0, chol0) 

# Clean the dataset

wcgs %>% mutate(bp_diff = sbp0 - dbp0) %>%
  select(id, sbp0, dbp0, bp_diff) %>%
  slice_min(bp_diff, n = 3)

wcgs %>% mutate(bp_diff = sbp0 - dbp0) %>%
  select(id, sbp0, dbp0, bp_diff) %>%
  filter(bp_diff >= 90)

# bp_diff: pulse pressure, typically in range 40-60. 
# These men did not have heart disease (in enrollment), high pulse pressure like this is unusual. 

wcgs <- wcgs %>%
  select(id, age0, bmi, sbp0, dbp0, smoker, behpat0, dibpat0, chol0, chd69, time169, arcus0, height0, weight0, ncigs0) %>%
  filter(complete.cases(.)) %>% 
  filter(bmi > 15) %>% 
  filter(sbp0 - dbp0 < 90) %>% 
  filter(chol0 <= 400)

View(wcgs)
# wcgs_ms: cleaned dataset

# Create additional variables
wcgs$behpat0f <- factor(wcgs$behpat0, levels = 1:4, label = c("A1", "A2", "B1", "B2"))
wcgs$dibpat0f <- factor(wcgs$dibpat0, levels = 0:1, label = c("B", "A"))
wcgs$smokerf <- factor(wcgs$smoker, levels = c(0, 1), labels = c("No", "Yes"))
wcgs$heightcm <- wcgs$height0 * 2.54
wcgs$weightkg <- wcgs$weight0 * 0.45359237

# wcgs$cholmmol <- wcgs$chol / 39
wcgs$time169y <- wcgs$time169 / 365.24
# Restrict follow-up to 5 years
wcgs$time169y5 <- pmin(wcgs$time169y, 5)
wcgs$chd695 <- (wcgs$chd69 == 1 & wcgs$time169y <= 5)

# View the structure of the dataset
str(wcgs)

# Create pair plot
pairs.panels(wcgs[, c("age0", "bmi", "weight0", "height0", "sbp0", "dbp0", "chol0", "chd69", "time169", "arcus0")])


# Explore distributions
ggplot(wcgs, aes(bmi)) + geom_histogram(position = "dodge", binwidth = 1) + labs(title = "BMI")
ggplot(wcgs, aes(chol0)) + geom_histogram(position = "dodge", binwidth = 1) + labs(title = "Cholesterol")

# Explore relationship between weight, smoking status, and CHD
# ggplot(wcgs, aes(x = weight0, y = ncigs0)) + geom_point(alpha = 0.2, position = position_jitter()) + facet_grid(~ chd69)
ggplot(wcgs, aes(x = weight0, y = chd69)) + geom_point(alpha = 0.2, position = position_jitter()) + facet_grid(~ smoker)
ggplot(wcgs, aes(x = bmi, y = chd69)) + geom_point(alpha = 0.2, position = position_jitter()) + facet_grid(~ arcus0)


# Test for linearity vs logit using boxTidwell
# p-value > 0.05 -> there is linearity -> linearity for all!

boxTidwell(formula = chd69 ~ age0, 
           other.x = ~ chol0 + sbp0 + dbp0 + bmi + smoker, 
           data = wcgs)

boxTidwell(formula = chd69 ~ chol0, 
           other.x = ~ age0 + sbp0 + dbp0 + bmi + smoker, 
           data = wcgs)

boxTidwell(formula = chd69 ~ sbp0, 
           other.x = ~ age0 + chol0 + dbp0 + bmi + smoker, 
           data = wcgs)

boxTidwell(formula = chd69 ~ bmi, 
           other.x = ~ age0 + chol0 + sbp0 + dbp0 + smoker, 
           data = wcgs)


# Build a logistic regression model (Model 1)
glm_model <- glm(chd69 ~ age0 + height0 + weight0 + chol0, data = wcgs, family = binomial)
summary(glm_model)


# Build another logistic regression model with BMI (Model 2)
glm_model_2 <- glm(chd69 ~ age0 + bmi + chol0, data = wcgs, family = binomial)
summary(glm_model_2)

# Run goodness-of-fit test for Model 2
with(summary(glm_model_2), 1 - deviance / null.deviance)


# Fit a logistic regression model with only BMI
lmod <- glm(chd69 ~ bmi, family = binomial, data = wcgs)
summary(lmod)


# https://www.statology.org/likelihood-ratio-test-in-r/?fbclid=IwAR0bJkF7DIBfG7NFv8VsC35NtDPua429a_wZ-saihqM1zS-rH2JUDkoMrpE
# LIKELIHOOD RATIO TEST
# null hypo: model lon hon va nho hon fit nhu nhau

# H0: b3 = b4 = b6 = 0
# H1: khong phai H0
# Nếu không reject H0: các b3, b4, b6 là insignificant -> bỏ được
# reject H0: b3, b4, b6 là significant -> không bỏ được ?

# Compare Model 2 and the BMI-only model using likelihood ratio test
lrtest(lmod, glm_model_2)

# Perform multicollinearity test (VIF) for Model 2
vif(glm_model_2)

# Drop diastolic blood pressure (dbp0) from Model 2 and perform likelihood ratio test
drop1(glm_model_2, test = "Chi")


# Perform logistic regression with BMI on the modified dataset (wcgs)
lmod <- glm(chd69 ~ bmi, family = binomial, data = wcgs)



# Add predicted probabilities and outcome predictions to the dataset

# Perform multiple logistic regression model with various predictors
glm_model_multi <- glm(chd69 ~ age0 + chol0 + sbp0 + dbp0 + bmi + smoker , data = wcgs, family = binomial)
summary(glm_model_multi)


vif(glm_model_multi)
# y = b0 + b1x1 + b2x2 

anova(glm_model_multi, test="Chi")
anova(glm_model_2, glm_model_multi, test = "Chi")
drop1(glm_model_multi, test="Chi")

# bỏ được dbp do p-value = 0.96 >>> 0.05
# tuy nhiên, chú ý: điều này không có nghĩa là dbp không ảnh hưởng tới chd

glm_model_3 <- glm(chd69 ~ dbp0, data=wcgs, family=binomial)
summary(glm_model_3)

# CORRELATION != CAUSATION
# Ta chi dang khao sat quan he giua cac bien, chứ không phải tính hệ quả do còn rất nhiều yếu tố ảnh hưởng tới bệnh mạch vành.
# model selection is not a stand-in for determining risk factors, let alone causative inference

# Stepwise backward model selection w.r.t AIC


# ban đầu có rất nhiều biến; với mỗi step ta vứt dần 1 biến có AIC thấp nhất cho đến khi AIC không thể giảm 

lmod <- glm(chd69 ~ age0 + height0 + weight0 + bmi + sbp0 + dbp0 + chol0 + dibpat0 + smoker + arcus0, family=binomial, wcgs)
summary(lmod)
lmodr <- step(lmod)

summary(lmodr)

vif(lmodr)
anova(lmodr, lmod, test="Chi")

drop1(lmodr, test="Chi")


# https://mattkcole.com/2017/01/22/the-problem-with-backward-selection/
# Nhận xét về phương pháp backward selection: 
# https://stats.stackexchange.com/questions/35353/why-applying-model-selection-using-aic-gives-me-non-significant-p-values-for-the 


# Tìm mô hình tốt nhất bằng duyệt trâu (exhaustive search)
wcgs_fit <- wcgs[, c("age0", "bmi", "height0", "weight0", "sbp0", "dbp0", "smoker", "chol0", "behpat0", "arcus0", "chd69")]
best <- bestglm(wcgs_fit, family=binomial, IC="BIC")
best

model1 <- glm(chd69 ~ age0 + sbp0 + smoker + chol0 + behpat0, data=wcgs_fit, family=binomial)
summary(model1)
vif(model1)


# Ta thử thêm các biến vào model tốt nhất 

# Reject H0 ở mức 5%
# Không reject H0 ở mức 1%

model2 <- glm(chd69 ~ age0 + sbp0 + smoker + chol0 + behpat0 + bmi, data=wcgs_fit, family=binomial)
summary(model2)
vif(model2)
lrtest(model2, model1)

model3 <- glm(chd69 ~ age0 + sbp0 + smoker + chol0 + behpat0 + height0 + weight0, data=wcgs_fit, family=binomial)
summary(model3)
vif(model3)
lrtest(model3, model1)

model4 <- glm(chd69 ~ age0 + sbp0 + smoker + chol0 + behpat0 + height0 + weight0 + arcus0, data=wcgs_fit, family=binomial)
summary(model4)
vif(model4)
lrtest(model4, model1)

model5 <- glm(chd69 ~ age0 + sbp0 + dbp0 + smoker + chol0 + behpat0 + height0 + weight0 + arcus0, data=wcgs_fit, family=binomial)
summary(model5)
vif(model5)
lrtest(model5, model1)


# Likelihood ratio test
lrtest(model5, model4, model3, model2, model1)
lrtest(model2)
lrtest(model1)

# (Chú ý, do model 1 so với các model khác đều reject null hypothesis ở mức 5% -> ta có thể ép về 1% để đảm bảo kết luận model 1 tốt)
# Nhận xét: ở độ tin cậy 1%, ta có thể kết luận model 1 tốt nhất trong 5 models ta khảo sát với likelihood ratio test.

# Survival analysis
table(wcgs$chd69)

survival <- Surv(wcgs$time169, wcgs$chd69)
head(survival)

# The first reading is 1664.0. The plus indicates that it was followed for 1664 days 
# without an event, so the measurement is impaired. However, the fifth reading
# 1885 has no plus indicating that the man developed heart disease after 1885 days.

kmfit <- survfit(Surv(time169,chd69)~1,data=wcgs)
km.plot <- ggsurvplot(kmfit,risk.table = T,break.time.by=10,xlim=c(0,110),
                      ylim=c(0.99,1)) 
km.plot 


kmfit <- survfit(Surv(time169,chd69)~1,data=wcgs)
km.plot <- ggsurvplot(kmfit,risk.table = T,break.time.by=30,xlim=c(0,365),
                      ylim=c(0.99,1)) 
km.plot 


kmfit.2 <- survfit(Surv(time169,chd69)~dibpat0f,data=wcgs)
km.plot.2 <- ggsurvplot(kmfit.2,risk.table = T,xscale=365.35,
                        break.time.by=365.25,ylim=c(0.8,1),tables.height=0.3) 
km.plot.2 

kmfit.2
summary(kmfit.2, times=c(365.25,5*365.25))
survdiff(Surv(time169,chd69)~dibpat0f,data=wcgs)

# ref: https://bggj.is/SurvivalAnalysis/introduction-to-survival-analysis.html

