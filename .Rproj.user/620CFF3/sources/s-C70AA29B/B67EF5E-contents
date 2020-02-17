### Libraries

library(caret)
library(fda)
library(fda.usc)
library(dplyr)

### Data
### Establish a clear set of inputs that can be used generally for statistical estimation.
### We can also start working to establish .Rdata files that have the data frame stored correctly.

data.raw <- read.table('lupustherm.txt')
temperatures <- seq(45, 90, 0.1)
thermogram.length <- length(temperatures)
thermograms.raw <- data.raw %>% select(1:thermogram.length)
class.raw <- data.raw %>% select(thermogram.length+1)


### thermogram.fda
### A function that takes thermogram data and produces a 'smoothed' version based on nbasis or lambda-pen
###
### x : a matrix (for thermograms this should be ...), accepta a data.frame also
### argvals : the temperatures
### basis.size : the n.basis used to produce fda representation
### 
### NEEDS: Derivative.

thermogram.fda <- function(x, argvals, nbasis)
{
  if(typeof(x)=='list') x<-as.matrix(x)
  thermogram.length <- length(argvals)
  samples.n <- dim(x)[1]
  
  fdata.raw <- fdata(x, argvals=argvals)
  reduced.basis.spline <- create.bspline.basis(range=range(argvals), nbasis = nbasis)
  tryCatch(Smoothing.Matrix <- S.basis(argvals, basis=reduced.basis.spline), 
           error=function(x) {Smoothing.Matrix<<-diag(thermogram.length)})
  
  ### Tryign to decide how to handle if the 'smoothing' is singular
  ### Need to just return an unsmoothed matrix, but its global and I can't find another way yet.
  
  fdata.smooth <- fdata.raw
  fdata.smooth$data <- fdata.raw$data%*%Smoothing.Matrix
  return(fdata.smooth)
}

fdata.working <- thermogram.fda(thermograms.raw, temperatures, 300)
# rm(Smoothing.Matrix)


### Fold Creator

fold.creator<-function(classes, folds=10, trials=10)
{
  folds.list<-list()
  for(j in  1:trials)
    folds.list[[j]]<-caret::createFolds(factor(classes), k=folds, list=FALSE)
  return(folds.list)
}

### Logistic Regression
### Function for performing LR and cross-validation.
### Returns the resulting predictions of test.x

lr.fit<-function(train.x, test.x, train.classes)
{
  glm.step<-glm(train.classes~., family="binomial", data=as.data.frame(train.x))
  return(glm.step)
}

model.chosen='lr.fit'

predictor.set <- fdata.working$data
classes <- class.raw$V452
folds <- 10
repeats <- 10
folds.list <- fold.creator(classes, folds, repeats)


if(model.chosen=='lr.fit') func<-lr.fit
func<-lr.fit
### loop over iterations
i=1

### loop over folds
fld=1


train.index <- which(folds.list[[i]]!=fld)

train.x <- predictor.set[train.index,]
test.x <- predictor.set[-train.index,]
train.classes <- classes[train.index]

model.result<-func(train.x, test.x, train.classes)

if(model.chosen=='lr.fit')
{
  coefficients <- model.result$coefficients[-1]
  predictions <- predict.glm(model.result, as.data.frame(test.x), type="response")
}

  



plot(predictions)

validator<-function(preds, classes, folds.list)
{
  
}

lr.kcv<-function(preds, classes, folds.list)
{
  sim.full<-preds
  output.predictions<-list()
  output.accuracy<-numeric()
  output.sensitivity<-numeric()
  output.specificity<-numeric()
  
  for(k in 1:trials)
  {
    trial.temp<-foreach(j=1:folds, .combine=append, .export="lr.fit") %dopar%
      {
        ### LENGTHS
        n.train<-length(which(folds.list[[k]]!=j))
        n.test<-length(which(folds.list[[k]]==j))
        
        ### CLASSES AND COUNTS
        train.classes<-classes[folds.list[[k]]!=j]
        test.classes<-classes[folds.list[[k]]==j]
        N0<-length(which(test.classes==0))
        N1<-length(which(test.classes==1))
        
        train.x<-sim.full[,folds.list[[k]]!=j]
        test.x<-sim.full[,folds.list[[k]]==j]
        
        n.train<-length(which(folds.list[[k]]!=j))
        n.test<-length(which(folds.list[[k]]==j))
        
        fit.out<-lr.fit(train.x, test.x, train.classes, test.classes)
        return(list(fit.out))
      }
    
    for(j in 1:folds)
    {
      output.predictions[[folds*(k-1)+j]]<-trial.temp[[j]]$predictions
      output.accuracy[folds*(k-1)+j]<-trial.temp[[j]]$accuracy
      output.sensitivity[folds*(k-1)+j]<-trial.temp[[j]]$sensitivity
      output.specificity[folds*(k-1)+j]<-trial.temp[[j]]$specificity 
    }
  }
  
  return(list(accuracy=output.accuracy, sensitivity=output.sensitivity, specificity=output.specificity, est.class.probs=output.predictions))
}