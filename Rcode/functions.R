library(gbm)
library(mgcv)
library(MASS)
library(rpart)
library(randomForest)


################################################################################
#-- Component Plots
################################################################################

#--- Things to Do
## add support for other classes: Date, POSIXct, ...
## Control by restricting edf, not sp. Also control categorical data.
## Add support for ordered factors (and convert numeric with few unique to this)

## still don't like the built in gam.mgcv smooths. They are not good
#  when x data has gaps
#   Note: tried using penalized b-splines with 50 knots and restricting min.sp=500
#         This min.sp=500 is rather arbitrary.
#   Want some way to specify maximum edf.



library(mgcv)

#-------------------------------------------------------
#-- function to plot all component plots
#-------------------------------------------------------
cp <- function(fmla,data,fam = poisson(),ylim=c(-2,2),max.df=6,plot=TRUE){ #,min.sp=100)
  data = model.frame(fmla,data)  # make correct transformations from fmla
  offset = as.vector(model.offset(data))
  y = model.response(data)    # response vector
  vars = attr(terms(data),"term.labels") # predictor variables
  p = length(vars)
  if(plot){
    op <- par(no.readonly=TRUE)
    on.exit(par(op))
    psqrt = sqrt(p)
    if(psqrt %% 1 == 0) xy = c(psqrt,psqrt)
    else xy = c(ceiling(psqrt),floor(psqrt))
    par(mfrow=c(xy[1],xy[2]),mar=c(2,2,1.5,1))
  }
  score = edf = numeric(p)
  for(j in 1:p){
    x = data[,vars[j]]
    fit = component.fit(x,y,offset,fam=fam,max.df=max.df,ylim=ylim,plot=plot)
    score[j] = fit$score
    edf[j] = fit$edf
    if(plot) title(paste0(vars[j]," (",round(score[j]),")"))
  }
return(data.frame(vars,score,edf,var.type=sapply(data[,vars],class),
                  stringsAsFactors=FALSE))  
}


#-------------------------------------------------------
#-- function to plot single smoothed fit
#-------------------------------------------------------
component.fit <- function(x,y,offset=NULL,fam=poisson(),max.df=6, #min.sp=500,
                          plot=TRUE, ylim=c(-2,2),xlab='',ylab='relative effect'){
## min.sp=500 sets edf<8 for all variables
## max.df=7 forces all numeric vars to have edf<6
# the points are log(y) - offset - intercept, which is similar to 
#   log(mu) = offset + intercept + smooth(x)
  
  
  if(plot){
    pt.color = "grey40"  #rgb(0,0,0,alpha=.8)
    color = col2rgb("grey80")/255
    CI.col = rgb(color[1],color[2],color[3],alpha=.6)
    line.col = "blue"
  }
  
  AIC.null = gam(y~1,family=fam,offset=offset)$aic   # AIC for baseline only model
  
  if(length(unique(x)) <= 3 & class(x) %in% c("numeric","integer")) {
    x = factor(x)
  }
  class.x = class(x)

  if(class.x %in% c("factor","character","logical")) {
    x = factor(x)
    gg = gam(y ~ s(x,bs = "re"),family = fam,offset = offset,method = "REML")
    
    if (plot) {
      plot(c(.5,nlevels(x) + .5),ylim,typ = 'n',xaxt = 'n',las = 1,
           xlab=xlab,ylab=ylab)
      axis(1,at = 1:nlevels(x))
      yy = fam$linkfun(y) - coef(gg)[1]
      if (!is.null(offset)) yy = yy - offset # note: log link
      points(jitter(as.numeric(x)),yy,pch = 19,cex = .4,col = pt.color)
      pp = predict(gg,data.frame(x = levels(x)),type = 'link',se.fit = TRUE)
      pp$fit = pp$fit - coef(gg)[1]  # remove intercept
      lower = pp$fit - 2 * pp$se.fit
      upper = pp$fit + 2 * pp$se.fit
      for (i in 1:nlevels(x)) {
        rect(i - .25,lower[i],i + .25,upper[i],col = CI.col)
        lines(i + c(-.25,.25),rep(pp$fit[i],2),col = line.col,lwd = 2)
      }
      abline(h = 0,lty = 3) #abline(h=log(mean(y)),lty=3)
      rug(jitter(as.numeric(x[yy < ylim[1]])),col = pt.color)
      rug(jitter(as.numeric(x[yy > ylim[2]])),side = 3,col = pt.color)
      box()
    }
  }

  if(class.x %in% c("numeric","integer")){
  nx = length(unique(x))

  #- Bspline with k=50 knots. 3rd degree spline, 2nd degree penalty (limits to liner fit)
  #  min.sp sets the minimum level of smoothing. Would rather set max df!
  #gg = gam(y~s(x,k=50,bs="ps",m=c(2,2)),family=fam,offset=offset,method="REML",
  #       min.sp=min.sp) 
  kmax = min(nx-1,max.df-1)
  gg = gam(y~s(x,k=kmax),family=fam,offset=offset,method="REML")  
  
  if(plot){
    yy = fam$linkfun(y) - coef(gg)[1]
    if(!is.null(offset)) yy = yy - offset # note: log link
    plot(x,yy,pch=19,cex=.4,ylim=ylim,las=1,typ='n',xlab=xlab,ylab=ylab)
    xseq = seq(min(x),max(x),length=100)
    pp = predict(gg,newdata=data.frame(x=xseq),type='link',se.fit=TRUE)
    pp$fit = pp$fit - coef(gg)[1]  # remove intercept
    lower = pp$fit - 2*pp$se.fit
    upper = pp$fit + 2*pp$se.fit
    polygon(c(xseq,rev(xseq)),c(lower,rev(upper)),col=CI.col,border=NA)
    points(x,yy,pch=19,cex=.4,col=pt.color)
    lines(xseq,pp$fit,col=line.col,lwd=2)    
    rug(jitter(x[yy<ylim[1]]),col = pt.color)         
    rug(jitter(x[yy>ylim[2]]),side=3,col = pt.color)  # show event that exceed bounds
    abline(h=0,lty=3) #abline(h=log(mean(y)),lty=3)   
    box()
    }
  }  
  score = AIC.null - gg$aic  # variable importance (componentwise)  
  edf = sum(gg$edf)
return(data.frame(score,edf))  
}


#-- This is old version:
component.plots <- function(fmla,data,plot=TRUE,fam = poisson(),ylim=c(-2,2),
                            max.df=6,min.sp=100, VARS="NONE",xy=0){
  
  data = model.frame(fmla,data)  # make correct transformations from fmla
  offset = as.vector(model.offset(data))
  y = model.response(data)    # response vector
  vars = attr(terms(data),"term.labels") # predictor variables
  p = length(vars)
  AIC.null = gam(y~1,family=fam,offset=offset)$aic   # AIC for baseline only model
  
  # calculate scores first and them order by score decreasingly "vars"
  AIC = numeric(p)
  for(j in 1:p){
    
    x = data[,vars[j]]
    if(length(unique(x)) <= 2 & class(x) %in% c("numeric","integer")){
      x = factor(x)
    }
    class.x = class(x)
    
    if(class.x %in% c("factor","character","logical")){
      x = factor(x)
      #gg = gam(y~x,family=fam,offset=offset)
      gg = gam(y~s(x,bs="re"),family=fam,offset=offset,method="REML")
    }
    
    if(class.x %in% c("numeric","integer")){
      nx = length(unique(x))
      #kmax = ifelse(nx < 10, nx-1, -1 ) # -1 lets gam choose kmax 
      kmax = min(nx-1,max.df-1)
      gg = gam(y~s(x,k=kmax),family=fam,offset=offset,method="REML")      
    }
    AIC[j] = gg$aic
  }
  score = AIC.null - AIC  # variable importance (componentwise)
  vars = vars[order(score, decreasing = T)]
  score = score[order(score, decreasing = T)]
  
  if(length(VARS) > 1) {
    VARS.index = which(vars %in% VARS)
    vars = vars[VARS.index]
    score = score[VARS.index]
    p = length(vars)
    panel.label = paste("(", letters[1:p], ")", sep = "")
  }
  
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  if(plot){ 
    if(length(xy) == 1) {
      psqrt = sqrt(p)
      if(psqrt %% 1 == 0) xy = c(psqrt,psqrt)
      else xy = c(ceiling(psqrt),floor(psqrt))
    }    
    par(mfrow=c(xy[1],xy[2]),mar=c(4,4,1.5,1))
    ylab = "" #"log(mu/offset)-Intercept" = smooth(x)
    color = col2rgb("grey80")/255
    CI.col = rgb(color[1],color[2],color[3],alpha=.6)
    line.col = "blue"
  }  
  
  for(j in 1:p){
    
    x = data[,vars[j]]
    if(length(unique(x)) <= 3 & class(x) %in% c("numeric","integer")){
      x = factor(x)
    }
    class.x = class(x)
    
    if(class.x %in% c("factor","character","logical")){
      x = factor(x)
      #gg = gam(y~x,family=fam,offset=offset)
      gg = gam(y~s(x,bs="re"),family=fam,offset=offset,method="REML")
      if(plot){
        plot(c(.5,nlevels(x)+.5),ylim,typ='n',xaxt='n',las=1,xlab='',ylab=ylab)
        axis(1,at=1:nlevels(x))
        yy = fam$linkfun(y) - coef(gg)[1]
        if(!is.null(offset)) yy = yy - offset # note: log link
        points(jitter(as.numeric(x)),yy,pch=19,cex=.4)       
        pp = predict(gg,data.frame(x=levels(x)),type='link',se.fit=TRUE)
        pp$fit = pp$fit - coef(gg)[1]  # remove intercept
        lower = pp$fit - 2*pp$se.fit
        upper = pp$fit + 2*pp$se.fit
        for(i in 1:nlevels(x)){
          rect(i-.25,lower[i],i+.25,upper[i],col=CI.col)
          lines(i+c(-.25,.25),rep(pp$fit[i],2),col=line.col,lwd=2)
        }
        lines(1:nlevels(x),pp$fit,col=line.col,lwd=2)          
        abline(h=0,lty=3) #abline(h=log(mean(y)),lty=3)  
        rug(jitter(as.numeric(x[yy<ylim[1]])))         
        rug(jitter(as.numeric(x[yy>ylim[2]])),side=3)        
        #bp = boxplot(y~x,col='grey80',ylim=c(0,25),las=1)
        #pp = predict(gg,data.frame(x=levels(x)),type='response') #,se.fit=TRUE)
        #lines(1:nlevels(x),pp,col=line.col,lwd=2)  
        #title(paste0(vars[j]," (",round(AIC.null-gg$aic),")"))
        #title(paste0(vars[j]," (",round(score[j]),")"))
        if(length(VARS) > 1) {
          mtext(side = 1, paste(panel.label[j], paste0(vars[j]," [",round(AIC.null-gg$aic),"]")), line = 2.5, cex = 1.2)
          mtext(side = 2, "Relative strength", line = 2, cex = 1.2)
        }
          
      }
    }
    
    if(class.x %in% c("numeric","integer")){
      nx = length(unique(x))
      
      #- Bspline with k=50 knots. 3rd degree spline, 2nd degree penalty (limits to liner fit)
      #  min.sp sets the minimum level of smoothing. Would rather set max df!
      min.smooth = 100  
      gg = gam(y~s(x,k=50,bs="ps",m=c(2,2)),family=fam,offset=offset,method="REML",
               min.sp=min.sp) 
      
      #kmax = ifelse(nx < 10, nx-1, -1 ) # -1 lets gam choose kmax 
      #kmax = min(nx-1,max.df-1)
      #gg = gam(y~s(x,k=kmax),family=fam,offset=offset,method="REML")      
    
      #kmax = max.df-1
      #gg = gam(y~s(x,k=kmax,bs="ps",m=c(2,1)),family=fam,offset=offset) 
            
      if(plot){
        #plot(gg,scheme=1,shift=coef(gg)[1],rug=FALSE,xlab="",ylim=ylim,las=1,
        #     shade.col=CI.col,col=line.col,lwd=2)       
        yy = fam$linkfun(y) - coef(gg)[1]
        if(!is.null(offset)) yy = yy - offset # note: log link
        plot(x,yy,pch=19,cex=.4,xlab="",ylab=ylab,ylim=ylim,las=1,typ='n')
        xseq = seq(min(x),max(x),length=100)
        pp = predict(gg,newdata=data.frame(x=xseq),type='link',se.fit=TRUE)
        pp$fit = pp$fit - coef(gg)[1]  # remove intercept
        lower = pp$fit - 2*pp$se.fit
        upper = pp$fit + 2*pp$se.fit
        polygon(c(xseq,rev(xseq)),c(lower,rev(upper)),col=CI.col,border=NA)
        points(x,yy,pch=19,cex=.4)
        lines(xseq,pp$fit,col=line.col,lwd=2)    
        rug(jitter(x[yy<ylim[1]]))         
        rug(jitter(x[yy>ylim[2]]),side=3)  # show event that exceed bounds
        abline(h=0,lty=3) #abline(h=log(mean(y)),lty=3)   
        # plot(gg,scheme=1,shift=coef(gg)[1],trans=fam$linkinv,ylim=log(range(y))) 
        # points(jitter(x),jitter(y))
        #title(paste0(vars[j]," (",round(AIC.null-gg$aic),")"))
        #title(paste0(vars[j]," (",round(score[j]),")"))
        if(length(VARS) > 1) {
          mtext(side = 1, paste(panel.label[j], paste0(vars[j]," [",round(score[j]),"]")), line = 2.5, cex = 1.2)
          mtext(side = 2, "Relative strength", line = 2, cex = 1.2)
        }          
      }
    }
  }
  return(data.frame(variable=vars,score=round(score)))
}




################################################################################
# Cross-Validation functions for intersection accident data
################################################################################

#== partition data for cross-validation
cvfolds <- function(n,k=10,seed) { 
  if(!missing(seed)) set.seed(seed)
  sample(rep(seq(k),length=n))
}



#== Cross-Val for Poisson Regression model
### Note: This isn't working with cross-validation. Perhaps due to some factors
#    whose values aren't included in one of the training sets. It seems that 
#    glm drops unused levels instead of giving a zero coefficient. 
#    Don't know what to do. Check out bestglm package?
#
#- TODO: Choose a reasonably short complexity/length to stop
cv.poisReg.old <- function(fmla,data,fold,show.pb=TRUE){
  X = model.matrix(fmla,data)
  Y = as.matrix(data[,as.character(fmla[[2]])])
  offset = model.offset(model.frame(fmla,data))  
  
  fit0 = glm.fit(X,Y,offset=offset,family=poisson())
  fit = glm(fmla,data=a,family=poisson)
  
  
  null.fmla = update(fmla,~1)  
  mu = numeric(nrow(data)) 
  K = sort(unique(fold))
  if(show.pb) pb = txtProgressBar(style=3,min=min(K),max=max(K))
  for(k in K) {
    test = which(fold == k)
    train = which(fold != k)
#fit0 = glm.fit(X[train,1],Y[train],offset=offset[train],family=poisson())
#fit = step(fit0,scope=fmla,direction="forward")
    
    
    fit0 = glm.my(fmla,data=data[train,],family=poisson)
    fit = step(fit0,trace=0)
    mu[test] = predict(fit,newdata=data[test,],type="response")
    if(show.pb) setTxtProgressBar(pb,k)
  }
  if(show.pb) close(pb)  
return(mu)  
}


cv.poisReg <- function(fmla,data,fold,max.complex=5,show.pb=TRUE){  
  X = model.matrix(fmla,data)
  Y = as.matrix(data[,as.character(fmla[[2]])])
  offset = model.offset(model.frame(fmla,data)) 
  
  tempdata = model.frame(fmla,data)  # make correct transformations from fmla
  vars = attr(terms(tempdata),"term.labels") # predictor variables

  
#   null.fmla = update(fmla,~1)  
  mu = matrix(NA,nrow(data),max.complex)
  K = sort(unique(fold))
  if(show.pb) pb = txtProgressBar(style=3,min=min(K),max=max(K))
  for(k in K) {
    test = which(fold == k)
    train = which(fold != k)
    
    ## Not complete
    fmla.rolling = as.formula(X5YrCrashCount~1) # add one variable to fmla.rolling each time
    fmla.scope = fmla # delete one variable from fmla.scope each time. NO NEED TO DO THIS?
    
    for(i in 1:max.complex) { # loop all the variables for now
      fit0 = glm.my(fmla.rolling,family=poisson,data=data[train,])  
      addProg = add1(fit0,scope=fmla)
      addVar = rownames(addProg[order(addProg$AIC), ][1, ])
      if(addVar == "<none>") next # if adding a variable doesn't help, just stop. No continue
      if (fmla.rolling[[3]] == 1) {
        fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), addVar), collapse = " "))  
      } else {
        fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), 
                                          fmla.rolling[[3]], "+", addVar), collapse = " "))  
      }
      if(!is.null(offset))
        fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), 
                                          fmla.rolling[[3]], "+ offset(log(Traffic))"), collapse = " ")) 
      
      #if(!is.null(offset)) fit.rolling = glm.my(fmla.rolling,data=data[train,],family=poisson, offset = offset[train])
      #else fit.rolling = glm.my(fmla.rolling,data=data[train,],family=poisson)
      fit.rolling = glm(fmla.rolling,data=data[train,],family=poisson) # bug when you mass around with "data"
      mu[test, i] = predict(fit.rolling,newdata=data[test,],type="response")
      #print(fmla.rolling)
      #print(i)
    }
    
    
    
    #fit1 = step(fit0,scope=fmla,direction="forward")
    
    #fit0 = glm.fit(X[train,1],Y[train],offset=offset[train],family=poisson())
    #fit = step(fit0,scope=fmla,direction="forward")
        
    #fit0 = glm.my(fmla,data=data[train,],family=poisson)
    #fit = step(fit0,trace=0)
    #mu[test] = predict(fit,newdata=data[test,],type="response")
    if(show.pb) setTxtProgressBar(pb,k)
  }
  if(show.pb) close(pb)  
  return(mu)  
} 

best.poisReg.fit <- function(fmla,data,cp.poisReg) {
  fmla.rolling = as.formula(X5YrCrashCount~1) # add one variable to fmla.rolling each time
  fmla.scope = fmla # delete one variable from fmla.scope each time. NO NEED TO DO THIS?
  
  for(i in 1:cp.poisReg) { # loop all the variables for now
    fit0 = glm.my(fmla.rolling,family=poisson,data=data)  
    addProg = add1(fit0,scope=fmla)
    addVar = rownames(addProg[order(addProg$AIC), ][1, ])
    if(addVar == "<none>") next # if adding a variable doesn't help, just stop. No continue
    if (fmla.rolling[[3]] == 1) {
      fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), addVar), collapse = " "))  
    } else {
      fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), 
                                        fmla.rolling[[3]], "+", addVar), collapse = " "))  
    }
    if(!is.null(offset))
      fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), 
                                        fmla.rolling[[3]], "+ offset(log(Traffic))"), collapse = " ")) 
    fit.rolling = glm(fmla.rolling,data=data,family=poisson) # bug when you mass around with "data"
  }
  return(fit.rolling)
}





#== Cross validation for negative bomonial
cv.negBino <- function(fmla,data,fold,max.complex=5,show.pb=TRUE){  
  X = model.matrix(fmla,data)
  Y = as.matrix(data[,as.character(fmla[[2]])])
  offset = model.offset(model.frame(fmla,data)) 
  
  tempdata = model.frame(fmla,data)  # make correct transformations from fmla
  vars = attr(terms(tempdata),"term.labels") # predictor variables
  
  
  #   null.fmla = update(fmla,~1)  
  mu = matrix(NA,nrow(data),max.complex)
  K = sort(unique(fold))
  if(show.pb) pb = txtProgressBar(style=3,min=min(K),max=max(K))
  for(k in K) {
    test = which(fold == k)
    train = which(fold != k)
    
    ## Not complete
    fmla.rolling = as.formula(X5YrCrashCount~1) # add one variable to fmla.rolling each time
    fmla.scope = fmla # delete one variable from fmla.scope each time. NO NEED TO DO THIS?
    
    for(i in 1:max.complex) { # loop all the variables for now
      fit0 = glm.nb.my(fmla.rolling,data=data[train,])  
      addProg = add1(fit0,scope=fmla)
      addVar = rownames(addProg[order(addProg$AIC), ][1, ])
      if(addVar == "<none>") next # if adding a variable doesn't help, just stop. No continue
      if (fmla.rolling[[3]] == 1) {
        fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), addVar), collapse = " "))  
      } else {
        fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), 
                                          fmla.rolling[[3]], "+", addVar), collapse = " "))  
      }
      if(!is.null(offset))
        fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), 
                                          fmla.rolling[[3]], "+ offset(log(Traffic))"), collapse = " ")) 
      
      fit.rolling = glm.nb.my(fmla.rolling,data=data[train,]) # bug when you mass around with "data"
      mu[test, i] = predict(fit.rolling,newdata=data[test,],type="response")
    }

    if(show.pb) setTxtProgressBar(pb,k)
  }
  if(show.pb) close(pb)  
  return(mu)  
}

best.negBino.fit <- function(fmla,data,cp.negBino) {
  fmla.rolling = as.formula(X5YrCrashCount~1) # add one variable to fmla.rolling each time
  fmla.scope = fmla # delete one variable from fmla.scope each time. NO NEED TO DO THIS?
  
  for(i in 1:cp.negBino) { # loop all the variables for now
    fit0 = glm.nb.my(fmla.rolling,data=data)  
    addProg = add1(fit0,scope=fmla)
    addVar = rownames(addProg[order(addProg$AIC), ][1, ])
    if(addVar == "<none>") next # if adding a variable doesn't help, just stop. No continue
    if (fmla.rolling[[3]] == 1) {
      fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), addVar), collapse = " "))  
    } else {
      fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), 
                                        fmla.rolling[[3]], "+", addVar), collapse = " "))  
    }
    if(!is.null(offset))
      fmla.rolling = as.formula(paste(c(fmla.rolling[[2]], as.character(fmla.rolling[[1]]), 
                                        fmla.rolling[[3]], "+ offset(log(Traffic))"), collapse = " ")) 
    
    fit.rolling = glm.nb.my(fmla.rolling,data=data) # bug when you mass around with "data"
  }
  return(fit.rolling)
}

cv.negBino.old <- function(fmla,data,fold,show.pb=TRUE){
  X = model.matrix(fmla,data)
  Y = as.matrix(data[,as.character(fmla[[2]])])
  offset = model.offset(model.frame(fmla,data))  
    
  null.fmla = update(fmla,~1)  
  mu = numeric(nrow(data)) 
  K = sort(unique(fold))
  if(show.pb) pb = txtProgressBar(style=3,min=min(K),max=max(K))
  for(k in K) {
    test = which(fold == k)
    train = which(fold != k)
    
    fit0 = glm.nb.my(fmla,data=data[train,])
    fit = step(fit0,trace=0)
    mu[test] = predict(fit,newdata=data[test,],type="response")
    if(show.pb) setTxtProgressBar(pb,k)
  }
  if(show.pb) close(pb)  
  return(mu)  
}



#== Cross-Val for glmnet model
cv.glmnetPois <- function(fmla,data,fold,alpha=0.8){
  X = model.matrix(fmla,data)[,-1]
  Y = as.matrix(data[,as.character(fmla[[2]])])
  offset = model.offset(model.frame(fmla,data))
  fit = cv.glmnet(X,Y,offset=offset,family="poisson",alpha=alpha,
                  foldid=fold,keep=TRUE)
  mu = exp(fit$fit.preval)
  ## This gives same result  
  #   lam.seq = glmnet(X,Y,offset=offset,family="poisson",alpha=alpha)$lambda
  #   mu = matrix(NA,nrow(data),length(lam.seq)) 
  #   K = sort(unique(fold))
  #   for(k in K) {
  #     test = which(fold == k)
  #     train = which(fold != k)  
  #     fit = glmnet(X[train,],Y[train],offset=offset[train],family="poisson",
  #                  alpha=alpha,lambda=lam.seq)    
  #     mu[test,] = predict(fit,newx=X[test,],offset=offset[test],type="response")
  #   }
  colnames(mu) = paste0("lambda_",1:ncol(mu))
  attr(mu,"lambda") = fit$lambda
  return(mu)  
}


#== Cross-Val for tree models
cv.treePois <- function(fmla, data, fold, cp.seq, show.pb = FALSE) { # other arguments use default of rpart()
  offset = model.offset(model.frame(fmla,data)) 
  if(missing(cp.seq)) { # default cp.seq is problematic. Not same for fmla and fmla.offset
    fit.full = rpart(fmla, data = data, method = "poisson", minsplit = 2,
                     minbucket = 2, cp = 0, xval = 0) # grow the full tree
    # cp.seq = seq(0, max(fit.full$cptable[, 1]), by = 0.003)
    cp.seq = fit.full$cptable[, 1]
  }
  mu = matrix(NA,nrow(data),length(cp.seq))
  K = sort(unique(fold))
  if(show.pb) pb = txtProgressBar(style=3,min=min(K),max=max(K))
  
  for(k in K) {
    test = which(fold == k)
    train = which(fold != k)
    fit.k = rpart(fmla, data = data[train,], method = "poisson", cp = 0,xval=0, minbucket=3)     
    for(i in 1:length(cp.seq)) {
      fit = prune(fit.k,cp=cp.seq[i])
      mu[test, i] =  predict(fit, newdata = data[test,], type = "vector")       
    }
    if(show.pb) setTxtProgressBar(pb,k)    
  }
  if(show.pb) close(pb)
  if(!is.null(offset)) mu = sweep(mu,1,exp(offset),"*")
  return(mu)
}


cv.treePois.old <- function(fmla,data,fold){
  offset = model.offset(model.frame(fmla,data))  
  mu = numeric(nrow(data))
  K = sort(unique(fold))
  for(k in K) {
    test = which(fold == k)
    train = which(fold != k)  
    fit = rpart(fmla,data=data[train,],method="poisson")
    mu[test] = predict(fit,newdata=data[test,],type="vector")
  } 
  if(!is.null(offset)) mu = mu * exp(offset)
  return(mu)
} 

#-- Cross-Val for boosting
cv.gbmPois <- function(fmla, data, fold, tree.seq=seq(500,10000,by=100),
                       interaction.depth = 3,..., show.pb=FALSE) {
  n.trees = max(tree.seq)
  mu = matrix(NA,nrow(data),length(tree.seq))
  K = sort(unique(fold))
  if(show.pb) pb = txtProgressBar(style=3,min=min(K),max=max(K))
  for(k in K) {
    test = which(fold == k)
    train = which(fold != k)
    fit = gbm(fmla, data=data[train,], distribution = "poisson", n.trees = n.trees,
               interaction.depth = interaction.depth,...) 
    mu[test,] = suppressWarnings(predict(fit, newdata=data[test,], n.trees = tree.seq, type = "response")) 
    if(show.pb) setTxtProgressBar(pb,k)
  }
  if(show.pb) close(pb)
  offset = model.offset(model.frame(fmla,data))  
  if(!is.null(offset)) mu = sweep(mu,1,exp(offset),"*")
  colnames(mu) = paste0("ntree_",tree.seq)
  return(mu)  
}

#-- Cross-Val for boosting (old version)
cv.gbmPois.old <- function(fmla, data, fold,max.trees=10000,
                       interaction.depth = 3,..., show.pb=FALSE) {
  mu = numeric(nrow(data))
  K = sort(unique(fold))
  if(show.pb) pb = txtProgressBar(style=3,min=min(K),max=max(K))
  for(k in K) {
    test = which(fold == k)
    train = which(fold != k)
    fit = gbm(fmla, data=data[train,], distribution = "poisson", n.trees = max.trees,
               interaction.depth = interaction.depth,...) 
    best.iter <- suppressWarnings(gbm.perf(fit, method="OOB",plot.it = FALSE))
    mu[test] = suppressWarnings(predict(fit, newdata=data[test,], n.trees = best.iter, type = "response"))
    if(show.pb) setTxtProgressBar(pb,k)
  }
  if(show.pb) close(pb)
  offset = model.offset(model.frame(fmla,data)) 
  if(!is.null(offset)) mu = mu * exp(offset)
  attr(mu,"best.iter") = best.iter
  return(mu)  
}




#== Returns Mean Absolute Error
mae <- function(mu,y){
  mu = as.matrix(mu)
  apply(mu,2,function(x) mean(abs(y-x)))
}

medae <- function(mu,y){
  mu = as.matrix(mu)
  apply(mu,2,function(x) median(abs(y-x)))
}

#-- Mean Squared Error
mse <- function(mu,y){
  mu = as.matrix(mu)
  apply(mu,2,function(x) mean((y-x)^2))
}

#-- Mean negative log likelihood (poisson distribution)
mlogL <- function(mu,y){
  mu = as.matrix(mu)
  apply(mu,2,function(x) -mean(dpois(y,x,log=TRUE)))
}



################################################################################
# Variable importance using AIC/likelihood and leave-one-out impact
# -Todo: chose an appropriate measure for all model
################################################################################
varImp.loo <- function(fmla, data, family, k = 5, seed=11032014) {
  data.temp = model.frame(fmla,data)  # make correct transformations from fmla
  offset = as.vector(model.offset(data.temp))
  vars = attr(terms(data.temp),"term.labels") # predictor variables
  p = length(vars)
  
  fold = cvfolds(nrow(data), k, seed) 
  ptm <- proc.time()
  
  score = numeric(p)
  if (model == "poisson") {
    cvError.null = suppressWarnings(cv.poisReg(fmla, data = a, fold = fold))
    
    for (j in 1:p) {
      if (!is.null(offset)) {
        fmla.temp = as.formula(paste("X5YrCrashCount ~ ", paste(paste(vars[-j], collapse = " + "), " + offset(log(Traffic))")))
      } else {
        fmla.temp = as.formula(paste("X5YrCrashCount ~ ", paste(vars[-j], collapse = " + ")))
      }       
      
      score[j] = glm(fmla.temp, data, family = "poisson")$aic - AIC.null
    }
  }
  else if (model == "treePois") {
      
  }
  proc.time() - ptm 
  
  
  dfscore = data.frame(variable=vars,score=round(score))
  dfscore = dfscore[order(dfscore$score, decreasing = T), ]
  return(dfscore)
}


##== Traditional variable importance
varImpStandard <- function(v, allVariables) { # scale variable "importance" from 0 to 100
  v <- abs(v) * 100 / max(abs(v)) # use absolute values when only coeffcients are provided
  
  f.temp <- function(x) { # strip the last digit for categorical variables
    if (substr(x, nchar(x), nchar(x)) %in% c(1:9)) return(substr(x, 1, nchar(x) - 1)) # strip 1 to 9?
    else return(x)
  }
  v.names.stripped <- unlist(lapply(names(v), FUN = f.temp))
  v.agg <- aggregate(v, by = list(v.names.stripped), FUN = max)
  v <- v.agg$x
  names(v) <- v.agg$Group.1 # collapse categorical dummy variables
  
  #if(!("log(TotalT5YearInM)" %in% names(v))) { # check whether offset is in it
  #  v <- c(v, 0)
  #  names(v) <- c(v.agg$Group.1, "log(TotalT5YearInM)1")
  #}
  
  v <- v[order(names(v))] # this is order by name. Can order by default order?
  v = v * 100 / sum(v) # standardized to add up to 100
  
  if(length(v) < length(allVariables)) {
    s <- rep(0, length(allVariables))
    names(s) <- allVariables
    matchedVariables <- allVariables[which(allVariables %in% names(v))]
    s[matchedVariables] <- v[matchedVariables]
    v = s
  } else if(length(v) > length(allVariables)) { # remove extra variables
    v = v[names(v) %in% allVariables]
  }
  
  return(v)
}

varImpStandard2 <- function(dfscore) {
  score = dfscore$score
  score = (score - min(score)) 
  score = score * 100 / max(score) 
  dfscore$score = score
  return(dfscore)
}

cp.rescale <- function(v, a, b) {# rescale v into (a, b)
  return((v - min(v)) * (b - a)/(max(v) - min(v)) + a)
}



