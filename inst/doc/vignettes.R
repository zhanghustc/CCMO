## ------------------------------------------------------------------------
library(CCMO)
data("SampleData",package="CCMO")
dim(SampleData)
head(SampleData)

## ------------------------------------------------------------------------
library(CCMO)
Y = SampleData[,1]
Gc = SampleData[,2]
Gm = SampleData[,12]
X = SampleData[,-(1:21)]
fit = singleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X)

## ------------------------------------------------------------------------
names(fit)

## ------------------------------------------------------------------------
fit$new

## ----eval=FALSE----------------------------------------------------------
#  fit$cov.new

## ----eval=FALSE----------------------------------------------------------
#  fit$log
#  fit$cov.log

## ------------------------------------------------------------------------
fit = OmnibusTest(fit,test=7:10)
fit$Omnibus

## ------------------------------------------------------------------------
Gc = SampleData[,2:11]
Gm = SampleData[,12:21]
system.time(fit1 <- multipleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X,cores=1))
system.time(fit2 <- multipleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X,cores=2))

## ----eval=FALSE----------------------------------------------------------
#  fit2[[2]]

## ------------------------------------------------------------------------
fit <- multipleSNP(Y,Gc,Gm,Xo=X,Xc=X,Xm=X,X.Gm=X,test=7:10)
fit[[2]]$Omnibus

