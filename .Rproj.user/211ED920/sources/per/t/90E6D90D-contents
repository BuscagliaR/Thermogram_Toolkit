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

samples.n <- dim(data.raw)[1]

thermograms.raw <- data.raw %>% select(1:thermogram.length)
class.raw <- data.raw %>% select(-(1:thermogram.length))

### Provide Names for any columns that are classifiers information.  So here, we indicate lupus/non-lupus

colnames(class.raw) <- 'SLE'
class.raw$SLE <- factor(class.raw$SLE)

### Visualization
### We can probably produce some nice visualization using ggplot.

# cols = rainbow(samples.n)
# plot(temperatures, thermograms.raw[1,], type='l', ylim = c(-0.25, 0.5), col=cols[1])
# for(i in 2:samples.n) lines(temperatures, thermograms.raw[i,], col=cols[i])

### Functional Data
### To be used for smoothing and derivatives

fdata.raw <- fdata(thermograms.raw, argvals=temperatures)
plot(fdata.raw)

######## Smoothing ############

### Choose a Basis

# full.spline <- create.bspline.basis(range=range(temperatures), nbasis = thermogram.length+2)

### We can tune from basis = 2 to length(argvals))
basis <- 450 ### Tuneable, if basis=length(arvals) then use raw.
reduced.basis.spline <- create.bspline.basis(range=range(temperatures), nbasis = basis+2)
Smoothing.Matrix <- S.basis(temperatures, basis=reduced.basis.spline)

fdata.smooth <- fdata.raw
fdata.smooth$data <- fdata.raw$data%*%Smoothing.Matrix
plot(fdata.smooth)

cbind(fdata.smooth$data[,1], fdata.raw$data[,1])

fdata.working <- fdata.smooth ### Working functional data set
toolkit.ws <- data.frame(Class = class.raw$SLE, X = fdata.working$data)  ### data.frame ws
colnames(toolkit.ws) <- c('Class', temperatures)

################################
### Data flow/cleaning is all above.
### Below is me starting to work on functions that should call the working data set.
################################


### A function that takes thermogram data and produces a 'smoothed' version based on nbasis and lambda-pen

### x : a matrix (for thermograms this should be ...)
### argvals : the temperatures
### basis.size : the n.basis used to produce fda representation


typeof(x)
typeof(thermograms.raw)

x <- thermograms.raw
argvals <- temperatures
nbasis <- length(argvals)-1

thermogram.fda <- function(x, argvals, nbasis)
{
  if(typeof(x)=='list') x<-as.matrix(x)
  thermogram.length <- length(argvals)
  samples.n <- dim(x)[1]
  
  fdata.raw <- fdata(x, argvals=argvals)
  reduced.basis.spline <- create.bspline.basis(range=range(argvals), nbasis = nbasis)
  Smoothing.Matrix <- S.basis(argvals, basis=reduced.basis.spline)
  
  fdata.smooth <- fdata.raw
  fdata.smooth$data <- fdata.raw$data%*%Smoothing.Matrix
  
  return(fdata.smooth)
}

fdata.temp <- thermogram.fda(thermograms.raw, temperatures, 5)
plot(fdata.temp)
for(j in seq(10, 200, 10)) plot(thermogram.fda(thermograms.raw, temperatures, j), main=j)
