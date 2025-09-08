setwd("~/Dropbox/Home/Desktop2/PhDProjects/FlowIs")
library(fields)
source("logisticRegression.R")
source("logisticRegressionSimple.R")
library( tictoc)
load("heart2020.rda")
cardio<- heart2020
cardio$HeartDisease<- as.integer(ifelse(cardio$HeartDisease =="Yes",1,0 ))
N<- nrow( cardio)

XDF<- cardio[,c("HeartDisease", 
                "BMI", "AgeCategory",
                "Smoking",         
               "AlcoholDrinking",  "Stroke","Sex")]
set.seed(222)
ind<- sample( 1:nrow(cardio), 3e5, replace=FALSE)
XDF<- XDF[ ind,] 
Age<- as.numeric( factor(XDF$AgeCategory ))
Age <- as.integer((Age-1)*5 +24 - 2)
X<- cbind( 1,
           XDF$BMI,
           Age,
           ifelse( XDF$Smoking =="Yes",1,0),
           ifelse( XDF$AlcoholDrinking =="Yes",1,0),
           ifelse( XDF$Stroke =="Yes",1,0),
           ifelse( XDF$Sex =="Male",1,0)
)

y<- XDF$HeartDisease

tic()
fit0<- glm( y~X-1, family=binomial())
cat( "number of fisher steps: ", fit0$iter, fill=TRUE)
toc()
tic()
fit1<- logisticRegression(X,y, tol=1e-7, verbose=FALSE)
cat( "number of fisher steps: ", fit1$iter, fill=TRUE)
toc()

tic()
fit2<- logisticRegressionSimple(X,y, tol=1e-7, verbose=FALSE)
cat( "number of fisher steps: ", fit1$iter, fill=TRUE)
toc()

test.for.zero(fit0$coefficients, fit1$beta, tag="betas", tol=1e-7)
test.for.zero(fit0$coefficients, fit2$beta, tag="betas", tol=1e-7)
tab<- summary(fit0)$coefficients
test.for.zero( tab[,2],fit1$SE, relative=FALSE, tag="SEs", tol=1e-5)
test.for.zero( tab[,2],fit2$SE, relative=FALSE, tag="SEs", tol=1e-5)


tic() 

fit1B<- logisticRegression(X,y, tol=1e-10, betaStart = fit1$beta*1.01,
                          verbose=FALSE)
cat( "number of fisher steps: ", fit1B$iter, fill=TRUE)
toc()
 

########## checking baseline option
indBase<- 2:3
baselineLogit <- X[,indBase]%*%fit1$beta[indBase]
betaSubset<- fit1$beta[-indBase]
XSubset<- X[,-indBase]
fit2<- logisticRegression(XSubset,y, betaStart= betaSubset*1.1,
                          baseline=baselineLogit,
                          tol=1e-10, verbose=FALSE)
test.for.zero(fit2$beta ,fit1$beta[-indBase], relative=FALSE, tag="Subset betas", tol=1e-7)


