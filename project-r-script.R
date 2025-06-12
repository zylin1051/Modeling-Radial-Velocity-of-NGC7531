###########################################
## Group Members: Z. LIN, H. LUO, J. ZHAO, W. ZHAO
###########################################



#################################################################################
library(caret)
library(plotly) 
library(mgcv)  
library(MASS)
library(knitr)
library(viridis)
library(ggplot2)  
options(warn=-1)
#################################################################################

## The galaxy data set can be found here :
##  https://hastie.su.domains/ElemStatLearn/
##  Go to the panel on the left -> Data -> Galaxy: (Info) (Data)

################## Import the Galaxy Data Set ###################################
galaxy <- read.csv("galaxy.txt")
font <- list(family = "Times New Roman", size = 10)
velocity <- galaxy$velocity
east.west <- galaxy$east.west
north.south <- galaxy$north.south
galaxy <- data.frame(y=galaxy$velocity, x1=galaxy$east.west, x2=galaxy$north.south)
n <- length(velocity)
p <- 2
y <- galaxy$y
galaxy.3d <- plot_ly(x = ~x1, y = ~x2, z = ~y, data=galaxy,
                     size = 20, type = "scatter3d", mode = "markers",
                     color = ~y)%>%layout(scene=list(xaxis = list(title="West-East"),
                                                     yaxis = list(title="South-North"),
                                                     zaxis = list(title="Velocity")),
                                          title = '', 
                                          font=font)
galaxy.3d <- galaxy.3d %>% colorbar(title = "Measured Radial Velocity")
galaxy.3d
#################################################################################



########################## Ordinary Least Squares ###############################
multilinear <- lm(y ~ x1 + x2, data = galaxy)

## Plot the fitted surface
multilinear.coef <- coef(multilinear)
ew <- seq(min(galaxy$x1), max(galaxy$x1), length.out = 100)
ns <- seq(min(galaxy$x2), max(galaxy$x2), length.out = 100)
multilinear.fit.vel <- outer(ew,ns,function(a,b){
  multilinear.coef[1]+multilinear.coef[2]*b+multilinear.coef[3]*a})
OLS.surface <- galaxy.3d %>% add_surface(x=~ew, y=~ns, z=~multilinear.fit.vel, 
                                         colors = "pink", opacity=0.5, showscale=FALSE,
                                         type = "surface", 
                                         mode = "markers", showlegend = F)
OLS.surface <- OLS.surface %>%layout(scene = list(xaxis = list(title="East-West"),
                                                  yaxis = list(title="North-South"),
                                                  zaxis = list(title="Velocity")),
                                     title = '') %>%hide_colorbar()
OLS.surface

## Compute AIC
multilinear.aic <- function(model, data){
  # model is linear model object
  # data is a data.frame containing the sample data
  y <- data$y
  X <- model.matrix(model)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  d <- sum(diag(H))
  ybar <- mean(y)
  yhat <- H%*%y
  sigmahat2 <- sum((y-yhat)^2)/(n-d)
  return(-2*mean(dnorm(y,H%*%y,sqrt(sigmahat2),log=TRUE)) + 2*mean(diag(H)))
}
multilinear.aic(multilinear, galaxy)
#################################################################################



############### Regularized Least Squares - L2 Norm Penalty #####################
## Compute AIC and the lambda that minimizes AIC
ridge.aic <- function(data, lambda){
  # data is a data.frame containing the sample data
  # lambda is the ridge penalty
  y <- data$y
  X <- cbind(1,data$x1,data$x2)
  p <- dim(X)[2]
  XtX <- crossprod(X)
  H <- X%*%solve(XtX + diag(lambda,p))%*%t(X)
  d <- sum(diag(H))
  ybar <- mean(y)
  yhat <- H%*%y
  sigmahat2 <- sum((y-yhat)^2)/(n-d)
  return(-2*mean(dnorm(y,H%*%y,sqrt(sigmahat2),log=TRUE)) + 2*mean(diag(H)))
}

## AIC vs. lambda
lambda.vals_ <- exp(seq(-5,5,length.out=100))
aic.vals_ <- sapply(lambda.vals_, function(k){ridge.aic(galaxy, k)})
plot(lambda.vals_, aic.vals_, main=expression(Ridge~Penalty:~AIC~vs.~lambda),
     xlab = expression(lambda), ylab = "AIC", type = "l", col = "red")
ridge.lambda.opt <- optimize(f=function(k){
  ridge.aic(galaxy,k)},interval=exp(seq(-5,5,length.out=100)))$minimum
ridge.aic.opt <- optimize(f=function(k){
  ridge.aic(galaxy,k)},interval=exp(seq(-5,5,length.out=100)))$objective
lines(rep(ridge.lambda.opt,2), c(0,20), lty = 2, col = "black")
ridge.lambda.opt
ridge.aic.opt

## Fit ridge regression model
galaxy.ridge <- function(lambda){
  y <- galaxy$y
  x1 <- galaxy$x1
  x2 <- galaxy$x2
  X <- cbind(1, x1, x2)
  XtX <- crossprod(X)
  p <- dim(X)[2]
  betahat <- solve(XtX + diag(lambda,p))%*%t(X)%*%y
  H <- X%*%solve(XtX + diag(lambda,p))%*%t(X)
  d <- sum(diag(H))
  yhat <- H%*%y
  sigmahat2 <- sum((y-yhat)^2)/(n-d)
  return(list(coefficients=betahat, sigmahat=sqrt(sigmahat2)))
}
ridge.coef <- galaxy.ridge(ridge.lambda.opt)$coefficients
ridge.fit.vel <- outer(ew,ns,function(a,b){
  ridge.coef[1]+ridge.coef[2]*b+ridge.coef[3]*a})
ridge.surface <- galaxy.3d %>% add_surface(x=~ew, y=~ns, z=~ridge.fit.vel, 
                                           colors = "magenta",opacity=0.5,
                                           type = "surface", mode = "markers", 
                                           showscale=FALSE, showlegend = F)
ridge.surface <- ridge.surface%>%layout(scene=list(xaxis=list(title="East-West"),
                                                   yaxis =list(title="North-South"),
                                                   zaxis = list(title="")),
                                          title ="") %>% hide_colorbar()
ridge.surface
#################################################################################



################## Polynomial Regression - Cubic polynomials ####################
poly.formula <- paste("y~x1+I(x1^2)+I(x1^3)",
                      "+x2+I(x2^2)+I(x2^3)",
                      sep="")
poly.formula <- as.formula(poly.formula)
poly <- lm(poly.formula, data=galaxy)
poly.coef <- coef(poly)
poly.fit.vel <- outer(ew,ns,function(a,b){
  poly.coef[1]+poly.coef[2]*b+poly.coef[3]*(b^2)+
    poly.coef[4]*(b^3)+poly.coef[5]*a+poly.coef[6]*(a^2)+poly.coef[7]*(a^3)
  })
poly.surface <- galaxy.3d %>% add_surface(x=~ew, y=~ns, z=~poly.fit.vel, 
                                          type = "surface", mode = "markers", 
                                          showlegend = F, colorscale = "Oranges",
                                          opacity=0.9, showscale=F)
poly.surface <- poly.surface%>%layout(scene = list(xaxis = list(title = "East-West"),
                                                   yaxis = list(title = "North-South"),
                                                   zaxis = list(title = "")),
                                      title = "")%>%hide_colorbar()
poly.surface

## Compute AIC
poly.aic <- function(model, data){
  # model is linear model object
  # data is a data.frame containing the sample data
  y <- data$y
  X <- model.matrix(model)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  d <- sum(diag(H))
  ybar <- mean(y)
  yhat <- H%*%y
  sigmahat2 <- sum((y-yhat)^2)/(n-d)
  return(-2*mean(dnorm(y,H%*%y,sqrt(sigmahat2),log=TRUE)) + 2*mean(diag(H)))
}
poly.aic(poly, galaxy)
#################################################################################




########################## Cubic Splines Additive Model #########################
## Choose lambda that minimizes AIC for the galaxy datasets
galaxy.spline.AIC <- function(data, lambda) {
  # data is a data.frame object containing the data
  # lambda is a scalar for the penalty lambda
  y <- data$y
  x1 <- data$x1
  x2 <- data$x2
  n <- length(y)
  ybar <- mean(y)
  # Point constraint is at the origin (0,0)
  add.smooth.1 <- smoothCon(s(x1, bs="bs", pc=0), data=data, absorb.cons = T)
  add.smooth.2 <- smoothCon(s(x2, bs="bs", pc=0), data=data, absorb.cons= T)
  smooth.X1 <- add.smooth.1[[1]]$X
  smooth.X2 <- add.smooth.2[[1]]$X
  smooth.S1 <- add.smooth.1[[1]]$S[[1]]
  smooth.S2 <- add.smooth.2[[1]]$S[[1]]
  smooth.S <- Matrix::bdiag(smooth.S1,smooth.S2)
  smooth.X <- cbind(smooth.X1, smooth.X2)
  smooth.XtX <- crossprod(smooth.X)
  smooth.Hat <- smooth.X %*% solve(smooth.XtX + lambda*smooth.S) %*% t(smooth.X)
  yhat <- smooth.Hat%*%(y-ybar) + ybar
  trHat <- sum(diag(smooth.Hat))
  sigmahat2 <- sum((y - yhat)^2)/(n-trHat)
  aic <- -2*mean(dnorm(y,yhat,sqrt(sigmahat2),log=TRUE)) + 2*mean(diag(smooth.Hat))
  return(aic) 
}

## AIC vs. lambda
lambda.vals <- exp(seq(-5,5,length.out=100))
aic.vals <- sapply(lambda.vals, function(k){galaxy.spline.AIC(galaxy, k)})
plot(lambda.vals, aic.vals, main=expression(Additive~Model:~AIC~vs.~lambda),
     xlab = expression(lambda), ylab = "AIC", type = "l", col = "red")
lambda.opt <- optimize(function(k){
  galaxy.spline.AIC(galaxy,k)}, interval = exp(seq(-5,5,length.out=100)))$minimum
aic.opt <- galaxy.spline.AIC(galaxy,lambda.opt) # AIC
lines(rep(lambda.opt,2), c(8,10), lty = 2, col = "black")
lambda.opt
aic.opt

## Fit model with lambda.opt and plot surface
galaxy.spline.pred <- function(lambda, newdata){
  # lambda is a scalar for the penalty lambda
  # newdata a data.frame containing the new covariates to do predictions
  y <- galaxy$y
  n <- length(y)
  ybar <- mean(galaxy$y)
  
  # Point constraint is chosen to be the sample mean of velocities
  add.smooth.1 <- smoothCon(s(x1, bs="bs", pc=0),data=galaxy,absorb.cons=T)
  add.smooth.2 <- smoothCon(s(x2, bs="bs", pc=0),data=galaxy,absorb.cons=T)
  smooth.X1 <- add.smooth.1[[1]]$X
  smooth.X2 <- add.smooth.2[[1]]$X
  smooth.S1 <- add.smooth.1[[1]]$S[[1]]
  smooth.S2 <- add.smooth.2[[1]]$S[[1]]
  smooth.S <- Matrix::bdiag(smooth.S1,smooth.S2)
  smooth.X <- cbind(smooth.X1, smooth.X2)
  smooth.XtX <- crossprod(smooth.X)
  smooth.betahat <- solve(smooth.XtX + lambda*smooth.S)%*%t(smooth.X)%*%(y-ybar)
  smooth.Hat <- smooth.X %*% solve(smooth.XtX + lambda*smooth.S) %*% t(smooth.X)
  yhat <- smooth.Hat%*%(y-ybar) + ybar
  trHat <- sum(diag(smooth.Hat))
  sigmahat2 <- sum((y - yhat)^2)/(n-trHat)
  
  # predictions
  new.x1 <- PredictMat(add.smooth.1[[1]],data.frame(x1 = newdata$x1))
  new.x2 <- PredictMat(add.smooth.2[[1]],data.frame(x2 = newdata$x2))
  new.X <- cbind(new.x1,new.x2)
  predictions <- new.X%*%smooth.betahat + ybar
  return(list(predictions=predictions,
              betahat=smooth.betahat,
              yhat=yhat,
              sigmahat2=sigmahat2))
}
new.coord <- expand.grid(x1=ew,x2=ns)
add.fit.vel <- galaxy.spline.pred(lambda.opt, new.coord)
add.fit.vel.pred <- matrix(add.fit.vel$predictions,ncol = 100,nrow = 100,byrow=T)
additive.surface <- galaxy.3d %>% add_surface(x=~ew, y=~ns, z=~add.fit.vel.pred, 
                                              type = "surface", mode = "markers", 
                                              showlegend = F, colorscale = 'Blues',
                                              opacity=0.9, showscale=F)
additive.surface <- additive.surface%>%layout(scene=list(xaxis=list(title="East-West"),
                                                         yaxis=list(title="North-South"),
                                                         zaxis=list(title="")),
                                              title = "")%>%hide_colorbar()
additive.surface
#################################################################################




############ Contrast with a Reference Point at the Origin (0,0) ################
galaxy.spline.contrast <- function(lambda, newdata, ref){
  # lambda is a scalar for the penalty lambda
  # newdata a data.frame containing the new covariates for predictions
  # ref is a data.frame containing the reference point
  y <- galaxy$y
  n <- length(y)
  ybar <- mean(galaxy$y)
  
  # Point constraint is chosen to be the sample mean of velocities
  add.smooth.1 <- smoothCon(s(x1, bs="bs", pc=0),data=galaxy,absorb.cons=T)
  add.smooth.2 <- smoothCon(s(x2, bs="bs", pc=0),data=galaxy,absorb.cons=T)
  smooth.X1 <- add.smooth.1[[1]]$X
  smooth.X2 <- add.smooth.2[[1]]$X
  smooth.S1 <- add.smooth.1[[1]]$S[[1]]
  smooth.S2 <- add.smooth.2[[1]]$S[[1]]
  smooth.S <- Matrix::bdiag(smooth.S1,smooth.S2)
  smooth.X <- cbind(smooth.X1, smooth.X2)
  smooth.XtX <- crossprod(smooth.X)
  smooth.betahat <- solve(smooth.XtX + lambda*smooth.S)%*%t(smooth.X)%*%(y-ybar)
  smooth.Hat <- smooth.X %*% solve(smooth.XtX + lambda*smooth.S) %*% t(smooth.X)
  yhat <- smooth.Hat%*%(y-ybar) + ybar
  
  # predictions
  new.x1 <- PredictMat(add.smooth.1[[1]],data.frame(x1 = newdata$x1))
  new.x2 <- PredictMat(add.smooth.2[[1]],data.frame(x2 = newdata$x2))
  ref.x1 <- PredictMat(add.smooth.1[[1]],data.frame(x1 = ref$x1))
  ref.x2 <- PredictMat(add.smooth.2[[1]],data.frame(x2 = ref$x2))
  new.X <- cbind(new.x1,new.x2)
  ref.X <- cbind(ref.x1, ref.x2)
  predictions <- new.X%*%smooth.betahat + ybar
  pred.ref <- ref.X%*%smooth.betahat + ybar
  contrasts <- predictions - rep(pred.ref, length(predictions))
  
  # 95% confidence intervals of contrasts
  smooth.Hat <- smooth.X %*% solve(smooth.XtX + lambda*smooth.S) %*% t(smooth.X)
  trHat <- sum(diag(smooth.Hat))
  sigmahat2 <- sum((y - yhat)^2)/(n-trHat)
  smooth.badConn <- solve(smooth.XtX + lambda*smooth.S)
  smooth.varbeta <- sigmahat2*(smooth.badConn%*%smooth.XtX%*%smooth.badConn)
  contrast.mat <- sweep(new.X,2,ref.X,'-')
  var.contrast <- contrast.mat%*%smooth.varbeta%*%t(contrast.mat)
  contrasts.se <- diag(sqrt(var.contrast))
  return(list(contrasts=contrasts,
              betahat=smooth.betahat,
              contrasts.se = contrasts.se))
}
origin <- data.frame(x1=0,x2=0)
add.cont.vel <- galaxy.spline.contrast(lambda.opt, new.coord, origin)
add.cont.vel.grid <- matrix(add.cont.vel$contrasts, ncol = 100, nrow = 100, byrow=T)
add.cont.vel.se <- matrix(add.cont.vel$contrasts.se, ncol = 100, nrow = 100, byrow=T)
contrast.surface <- plot_ly(x=~ew, y=~ns, z=~add.cont.vel.grid, type = "surface", 
                            mode = "markers", showlegend = F, color=~add.cont.vel.grid,
                            opacity=0.9, 
                            showscale=T)%>%colorbar(title="",len=0.75, tickness=20)
contrast.surface%>%layout(scene = list(xaxis = list(title = "West-East"),
                                       yaxis = list(title = "South-North"),
                                       zaxis = list(title = "")),title = "")
contrast.surface %>% hide_colorbar()
image(ew, ns, -add.cont.vel.grid,
      col = viridis(100),
      xlab = "East-West",
      ylab = "South-North",
      main = "Heat Map of Contrast")
contour(x = ew, y = ns, z = -add.cont.vel.grid, nlevels = 50, xlab = "East-West",
        ylab = "South-North", main = "", add = TRUE)
lines(c(-25, 25), c(45, -45), lty = 2, lwd = 2, col = "red")
text(8.5, -19, 
     labels = "Steepest Change", col = "red", srt=-30, fonts="Arial")
#################################################################################



############################## Residuals Diagnostics ############################
result <- galaxy.spline.pred(lambda.opt,galaxy)
yhat <- result$predictions
sigmahat2 <- result$sigmahat2
residuals <- (y - yhat)/sqrt(sigmahat2)

# Visualize residuals
ggplot(data.frame(x1 = galaxy$x1, x2 = galaxy$x2, residuals = residuals),
       aes(x = x1, y = residuals)) + 
  ggtitle("Studentized Residuals vs Fitted Values") +
  geom_point() +
  xlab("Fitted Values") +
  ylab("Residuals") +
  geom_hline(yintercept = c(-1.96, 0, 1.96), color = "red", lty = 2)  

qqnorm(residuals, main = "Residuals Normal Q-Q Plot")
qqline(residuals, col = "red") 

hist(residuals, breaks = 20, probability = T, main = "Histogram of Residuals",
     xlab = "Studentized Residuals")
rr <- seq(-5, 5, length.out = 1000)
dd <- dnorm(rr)
lines(rr, dd, col = "red", lty = 1)
legend("topleft", "N(0,1) PDF", lty = 2, col = "red")
#################################################################################



########################## Cross Validation #####################################
## Function for cross-validation for the additive model
additive.cv <- function(train, test, lambda=lambda.opt){
  # train and test are data.frame objects containing the train and test set
  y.train <- train$y
  y.test <- test$y
  ybar.train <- mean(train$y)
  
  # Point constraint chosen to be the mean of the train data
  add.smooth.1 <- smoothCon(s(x1, bs="bs", pc=0),data=train,absorb.cons=T)
  add.smooth.2 <- smoothCon(s(x2, bs="bs", pc=0),data=train,absorb.cons=T)
  smooth.X1 <- add.smooth.1[[1]]$X
  smooth.X2 <- add.smooth.2[[1]]$X
  smooth.S1 <- add.smooth.1[[1]]$S[[1]]
  smooth.S2 <- add.smooth.2[[1]]$S[[1]]
  smooth.S <- Matrix::bdiag(smooth.S1,smooth.S2)
  smooth.X <- cbind(smooth.X1, smooth.X2)
  smooth.XtX <- crossprod(smooth.X)
  smooth.betahat <- solve(smooth.XtX+lambda*smooth.S)%*%t(smooth.X)%*%(y.train-ybar.train)
  
  # predictions
  new.x1 <- PredictMat(add.smooth.1[[1]],data.frame(x1 = test$x1))
  new.x2 <- PredictMat(add.smooth.2[[1]],data.frame(x2 = test$x2))
  new.X <- cbind(new.x1,new.x2)
  predictions <- new.X%*%smooth.betahat + ybar.train
  return(predictions)
}

## Function that partitions the data
partition.Kfold <- function(k, data, seed=4446){
  # k is the number of folds
  # data is a data.frame object containing the data
  # seed is a number for the random seed
  set.seed(seed)
  n <- dim(data)[1]
  if (k==n){ # leave-one-out if k==n
    partitions <- list()
    for (i in 1:n){
      partitions[[i]] <- data[i,]
    }
    return(partitions)
  }
  i <- c(1:n)
  size <- floor(n/k)
  partitions <- list()
  for (j in 1:k){
    if (j==k){# if k does not divide n, the last fold will have a bit more data
      partitions[[j]] <- data[i,]
    }else{
      selects <- sample(i, size, replace=F)
      remove <- sapply(selects, function(k){which(i==k)})
      i <- i[-remove]
      chosens <- data[selects, ]
      partitions[[j]] <- chosens
    }
  }
  return(partitions)
}

## cross-validation
cross.validate.project <- function(data.par, mod='additive') {
  # data.par is an object returned by partition.Kfold
  # mod is a type of models used in this project: "lm","poly","ridge","additive"
  stopifnot(mod%in%c('lm','poly','ridge','additive'))
  n <- sum(unlist(lapply(data.par, function(k){dim(k)[1]})))
  k <- length(data.par)
  loss <- 0
  for (i in 1:k){
    test <- data.par[[i]]
    train <- Reduce(rbind, data.par[-i])
    y.test <- test$y
    ## ordinary linear regression
    if (mod=="lm"){
      model <- lm(y~x1+x2, data=train)
      yhat.test <- predict(model, test)
    } 
    ## polynomial regression
    else if (mod=="poly"){
      poly.formula <- paste("y~x1+I(x1^2)+I(x1^3)",
                            "+x2+I(x2^2)+I(x2^3)",sep="")
      poly.formula <- as.formula(poly.formula)
      model <- lm(poly.formula, data=train)
      yhat.test <- predict(model, test)
    }
    else if (mod=="ridge"){
      y.train <- train$y
      x1.train <- train$x1
      x2.train <- train$x2
      x1.test <- test$x1
      x2.test <- test$x2
      X.train <- cbind(1, x1.train, x2.train)
      X.test <- cbind(1, x1.test, x2.test)
      XtX.train <- crossprod(X.train)
      p.train <- dim(X.train)[2]
      betahat.train <- solve(XtX.train + diag(ridge.lambda.opt,p.train))
      betahat.train <- betahat.train%*%t(X.train)%*%y.train
      yhat.test <- X.test%*%betahat.train
    } else{
      yhat.test <- additive.cv(train, test)
    }
    loss.i <- sum((y.test-yhat.test)^2)
    loss <- loss + loss.i
  }
  return(sum(loss)/n)
}

## Perform cross validations on the models
folds <- 10
galaxy.par <- partition.Kfold(folds, galaxy)
galaxy.par.loo <- partition.Kfold(n, galaxy)
CV.kfold.score <- c(cross.validate.project(galaxy.par,"lm"),
                    cross.validate.project(galaxy.par,"ridge"),
                    cross.validate.project(galaxy.par,"poly"),
                    cross.validate.project(galaxy.par,"additive"))
CV.kfold.score <- sapply(CV.kfold.score, function(s){round(s, 3)})
CV.loo.score <- c(cross.validate.project(galaxy.par.loo,"lm"),
                  cross.validate.project(galaxy.par.loo,"ridge"),
                  cross.validate.project(galaxy.par.loo,"poly"),
                  cross.validate.project(galaxy.par.loo,"additive"))
CV.loo.score <- sapply(CV.loo.score, function(s){round(s, 3)})
#################################################################################



########################### Model Performance Summary ###########################
train.err.OLS <- mean((galaxy$y-predict(multilinear, galaxy))^2)
train.err.ridge <- mean((galaxy$y-(ridge.coef[1]+ridge.coef[2]*galaxy$x1+ridge.coef[3]*galaxy$x2))^2)
train.err.poly <- mean((galaxy$y-predict(poly, galaxy))^2)
train.err.add <- mean((galaxy$y-galaxy.spline.pred(lambda.opt, 
                                                   data.frame(x1=0,x2=0))$yhat)^2)
train.err <- round(rbind(train.err.OLS, train.err.ridge, train.err.poly, train.err.add),3)
AIC.score <- c(multilinear.aic(multilinear, galaxy),
               ridge.aic.opt,
               poly.aic(poly, galaxy),
               aic.opt)
AIC.score <- sapply(AIC.score, function(s){round(s, 3)})
rows <- c("OLS", "Ridge Penalty", 
          "Cubic Polynomial",
          "Additive")
summary <- cbind(rows, AIC.score,  train.err, CV.kfold.score, CV.loo.score)
columns <- c("Method", "AIC", "Training Error",
             "Est. MPSE: K(=10)-Fold CV", "Est. MPSE: leave-one-out CV")
sum.table <- kable(summary, col.names = columns, booktab=T, align="l")
sum.table 
#################################################################################




############################## Plots and Tables #################################

## Overview of the Data
galaxy.3d  


## OLS Fitted Surface
OLS.surface 


## Ridge regression
# fitted surface
ridge.surface 
# AIC vs. Lambda
plot(lambda.vals_, aic.vals_, main=expression(Ridge~Penalty:~AIC~vs.~lambda),
     xlab = expression(lambda), ylab = "AIC", type = "l", col = "red")
lines(rep(ridge.lambda.opt,2), c(0,20), lty = 2, col = "black")


## Polynomial Regression Fitted Surface
poly.surface


## Additive Model 
# fitted surface
additive.surface 
# AIC vs. Lambda
plot(lambda.vals, aic.vals, main=expression(Additive~Model:~AIC~vs.~lambda),
     xlab = expression(lambda), ylab = "AIC", type = "l", col = "red")
lines(rep(lambda.opt,2), c(8,10), lty = 2, col = "black")


## Contrast surface
# with colourbar
contrast.surface%>%layout(scene = list(xaxis = list(title = "West-East"),
                                       yaxis = list(title = "South-North"),
                                       zaxis = list(title = "")),title = "")
# without colourbar
contrast.surface %>% hide_colorbar()


## Contrast Heatmap
image(ew, ns, -add.cont.vel.grid,
      col = viridis(100),
      xlab = "East-West",
      ylab = "South-North",
      main = "Heat Map of Contrast")
contour(x = ew, y = ns, z = -add.cont.vel.grid, nlevels = 50, xlab = "East-West",
        ylab = "South-North", main = "", add = TRUE)
lines(c(-25, 25), c(45, -45), lty = 2, lwd = 2, col = "red")
text(8.5, -19, 
     labels = "Steepest Change", col = "red", srt=-30, fonts="Arial")


## Residuals versus Fitted Values
ggplot(data.frame(x1 = galaxy$x1, x2 = galaxy$x2, residuals = residuals),
       aes(x = x1, y = residuals)) + 
  ggtitle("Studentized Residuals vs Fitted Values") +
  geom_point() +
  xlab("Fitted Values") +
  ylab("Residuals") +
  geom_hline(yintercept = c(-1.96, 0, 1.96), color = "red", lty = 2)  


## Residuals Q-Q Plot
qqnorm(residuals, main = "Residuals Normal Q-Q Plot")
qqline(residuals, col = "red") 


## Histogram of Residuals
hist(residuals, breaks = 20, probability = T, main = "Histogram of Residuals",
     xlab = "Studentized Residuals")
rr <- seq(-5, 5, length.out = 1000)
dd <- dnorm(rr)
lines(rr, dd, col = "red", lty = 1)
legend("topleft", "N(0,1) PDF", lty = 2, col = "red")


## Summary Table for Models
sum.table 


############################# END OF SCRIPT #####################################

