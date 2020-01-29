######### Function #############

## Classification using penalized LDA

library('penalizedLDA')
library("pROC")


##### Data Import #####

### This data is created using Data_Import_Cleaning where we can choose smoothing and derivatives.
plot(fdata.working) ### Working functional data set


##### Setup #####

x <- fdata.working$data
y <- as.numeric(class.raw$SLE)
set.seed(10)
folds.all <- caret::createFolds(class.raw$SLE)  ### prefer as input or we randomly generate the folds
lambda.seq <- exp(seq(-6, 2, 0.01))  ### This could be an input, best to suggest evalauating this outide of functions.

##### Training #####

fold.id <- 2

folds.all[[fold.id]]

x.test <- x[folds.all[[fold.id]],]
y.test <- y[folds.all[[fold.id]]]

x.train <- x[-folds.all[[fold.id]],]
y.train <- y[-folds.all[[fold.id]]]

penlda.cv.out <- PenalizedLDA.cv(x=x.train, y=y.train, lambdas=lambda.seq, nfold = 10)
### plot(penlda.cv.out)

lambda.set <- penlda.cv.out$bestlambda ### Record the lambda used
penlda.out <- PenalizedLDA(x=x.train, y=y.train, lambda=lambda.set, K=1)

discrim.temp <- penlda.out$discrim
### plot(seq(45, 90, 0.1), penlda.out$discrim)  ### Want to save this vector

##### Test Set #####

xte <- scale(x.test, center = apply(x.train, 2, mean), scale = penlda.out$wcsd.x) #normalize test data
probs.log.temp <- xte%*%penlda.out$discrim #calculte posterior log-odds
probs.temp <- 1/(1+exp(-probs.log.temp))



roc.temp <- pROC::roc(y.test, as.vector(probs.temp))
# plot(my.roc)
auc.temp <- as.numeric(my.roc$auc)

##### END #####



















penlda.cv.out <- PenalizedLDA.cv(x=x, y=y)


penlda.cv.out <- PenalizedLDA.cv(x=x, y=y, lambdas=lambda.seq, nfold = 10)

plot(penlda.cv.out)
penlda.cv.out$bestlambda

### The start of a function

set.seed(10)



### I have figured out how to get posterior estimates of the probability for the binary cases.

fold.id <- 2

folds.all[[fold.id]]

x.test <- x[folds.all[[fold.id]],]
y.test <- y[folds.all[[fold.id]]]

x.train <- x[-folds.all[[fold.id]],]
y.train <- y[-folds.all[[fold.id]]]

penlda.cv.out <- PenalizedLDA.cv(x=x.train, y=y.train, lambdas=lambda.seq, nfold = 10)

lambda.set <- penlda.cv.out$bestlambda

penlda.out <- PenalizedLDA(x=x.train, y=y.train, lambda=lambda.set, K=1)

plot(seq(45, 90, 0.1), penlda.out$discrim)  ### Want to save this vector

# penlda.out$wcsd.x # group deviations for normalizing predictors

xte <- scale(x.test, center = apply(x.train, 2, mean), scale = penlda.out$wcsd.x) #normalize test data
probs <- xte%*%penlda.out$discrim #calculte posterior
probs
predict.test <- probs>0
factor(predict.test, c(TRUE, FALSE), c(0,1))
y.test-1

table(y.test-1, predict.test)
mean(y.test-1==predict.test)

cbind(probs, penlda.preds$ypred)

penlda.preds <- predict(penlda.out, x.test)
acc <- mean(penlda.preds$ypred==y.test)  ### There does not seem to be a way to get posterior probabilities, so we may not include this
acc

predict.penlda

scale(x.test, mean(x.train), sd(x.train))

### A loop

set.seed(10)

folds.n <- 10
folds.all <- caret::createFolds(class.raw$SLE, k = folds.n)
lambda.seq <- exp(seq(-6, 2, 0.01))  ### This could be an input, best to suggest evalauating this outide of functions.
acc <- vector()

for(fold.id in 1:folds.n)
{
  x.test <- x[folds.all[[fold.id]],]
  y.test <- y[folds.all[[fold.id]]]
  
  x.train <- x[-folds.all[[fold.id]],]
  y.train <- y[-folds.all[[fold.id]]]
  
  penlda.cv.out <- PenalizedLDA.cv(x=x.train, y=y.train, lambdas=lambda.seq, nfold = 10)
  
  lambda.set <- penlda.cv.out$bestlambda
  
  penlda.out <- PenalizedLDA(x=x.train, y=y.train, lambda=lambda.set, K=1)
  
  penlda.preds <- predict(penlda.out, x.test)
  acc[fold.id] <- mean(penlda.preds$ypred==y.test)  ### There does not seem to be a way to get posterior probabilities, so we may not include this
}

boxplot(acc)  ### Got about 72% mean test set accuracy in first run, on a smoothed set
## Second run was on original data, seemed different, tuning should be interesting.

### Enforcing a fused penalty
### This is something of interest but takes double the tuning (so probably quadruple the effort)

# set.seed(10)
# 
# folds.n <- 10
# folds.all <- caret::createFolds(class.raw$SLE, k = folds.n)
# lambda.seq <- exp(seq(-6, 2, 0.01))  ### This could be an input, best to suggest evalauating this outide of functions.
# acc <- vector()
# 
# for(fold.id in 1:folds.n)
# {
#   x.test <- x[folds.all[[fold.id]],]
#   y.test <- y[folds.all[[fold.id]]]
#   
#   x.train <- x[-folds.all[[fold.id]],]
#   y.train <- y[-folds.all[[fold.id]]]
#   
#   penlda.cv.out <- PenalizedLDA.cv(x=x.train, y=y.train, lambdas=lambda.seq, nfold = 10, type='ordered')
#   
#   lambda.set <- penlda.cv.out$bestlambda
#   
#   penlda.out <- PenalizedLDA(x=x.train, y=y.train, lambda=lambda.set, K=1, type='ordered')
#   
#   penlda.preds <- predict(penlda.out, x.test)
#   acc[fold.id] <- mean(penlda.preds$ypred==y.test)  ### There does not seem to be a way to get posterior probabilities, so we may not include this
# }
# 
# boxplot(acc)  ### Got about 72% mean test set accuracy in first run
