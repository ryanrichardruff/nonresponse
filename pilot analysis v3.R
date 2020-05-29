#Microbiome nonresponder analysis
#Clear environment
rm(list = ls())

packages <- c('coefplot', 'plotmo', 'glmnet')
install.packages(packages)

library(glmnet)
library(plotmo)
library(coefplot)

options(max.print=1000000)

library(readxl)
otu <- read_excel("SDF OTU full data table for analysis RUFF _ analytic set.xlsx")
View(otu)   



#isolate amongst saliva, then create two subesets. Run program twice with the different subsetted "newdata" for both models
#current two comparisons: nonresponders to everyone else; nonresponders to SDF
newdata <- subset(otu, site==1)
#newdata <- subset(otu, site==1 & (group=="non_responder" | group=="sdf"))

#isolate all OTUs from original set
firstcol= which(colnames(newdata)=="NewCleanUpReferenceOTU3841")
lastcol= which(colnames(newdata)=="Vdisp")

#remove shared column names
names(newdata) <- sub("NewCleanUpReference","", names(newdata))
names(newdata) <- sub("NewReference","", names(newdata))
names(newdata)

#remove OTUs with at least 75% as 0 for all patient saliva sample, restricts variable list. Currently not using this
smallset<-newdata[, c(firstcol:lastcol)]
smallset <- smallset[colMeans(smallset==0) <= 0.75]


#cv glmnet
xmat<-newdata[, c(firstcol:lastcol)]
y_var1 <- newdata$nonres
x_var1 <- data.matrix(xmat)
lseq <- 10^seq(-3, 5, length.out=100)
lasso_cv <- cv.glmnet(x_var1, y=as.factor(y_var1), alpha=1, lambda=lseq, standardize=TRUE, nfolds=10, family="binomial")
plot(lasso_cv)
#save output to file
sink('analysis-outputv2.txt')
coef(lasso_cv)
sink() #close save to output

coefplot(lasso_cv)
coefpath(lasso_cv, showLegend="onmouseover")


#regular lasso based on cv glmnet
lcv <-lasso_cv$lambda.min
fit4<-glmnet(x_var1, y=as.factor(y_var1), lambda=lcv, family="binomial")
sink('analysis-lambda minv2.txt')
predict(fit4,type="response",newx=x_var1[2:5,])
predict(fit4,type="coefficients")
sink()
plot(fit4, xvar="lambda", label=TRUE)
plot_glmnet(fit4, xvar="norm", label=TRUE, s=lcv)

coefplot(fit4, lambda=lcv, sort='magnitude')
coefpath(fit4)

























#glnet with forcing lambda to be what stata reports, keeps same vars but not same numbers (algorithm?)
lcv <-lasso_cv$lambda.min
model_cv <- glmnet(x_var1, y=as.factor(y_var1), alpha=1, lambda=.20739863, standardize=TRUE, family="binomial")
coef(model_cv)
plot(model_cv)








#lasso
#https://www.rstatisticsblog.com/data-science-in-action/lasso-regression/
y_var <- newdata$nonres
x_var <- data.matrix(smallset)
#x_var <- data.matrix(smallset[, c(fcol:lcol)])  #dont need this since smallset only includes predictor vars

#optimal lambda
lambda_seq <- 10^seq(2, -2, by = -.1)
cv_output <- cv.glmnet(x_var, y=as.factor(y_var), alpha = 1, lambda = lambda_seq, family="binomial")

# identifying best lamda
best_lam <- cv_output$lambda.min

# Rebuilding the model with best lamda value identified
lasso_best <- glmnet(x_var, y=as.factor(y_var), alpha = 1, lambda = best_lam, family="binomial")
pred <- predict(lasso_best, s = best_lam, newx = x_var)

#Get retained variables
coef(lasso_best)

# Plot variable coefficients vs. shrinkage parameter lambda.
plot(lasso_best, pch=19)












