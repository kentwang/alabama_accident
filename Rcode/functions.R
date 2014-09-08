## Include offsets,
## Order plots by score
## the ylim for class=numeric is fixed. Add to function input.
## add support for other classes: Date, POSIXct, ...


library(mgcv)


component.plots <- function(fmla,data,plot=TRUE,fam = poisson(),yrng=c(0,25)){
  # data = model.frame(fmla,data)
  vars = all.vars(fmla)
  y = data[,vars[1]]
  vars = vars[-1]
  p = length(vars)
  AIC.null = gam(y~1,family=fam)$aic   # AIC for baseline only model
  
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  if(plot){ 
    psqrt = sqrt(p)
    if(psqrt %% 1 == 0) xy = c(psqrt,psqrt)
    else xy = c(ceiling(psqrt),floor(psqrt))
    par(mfrow=c(xy[1],xy[2]),mar=c(2,4,1.5,1))
  }  
  
  AIC = numeric(p)
  for(j in 1:p){
    x = data[,vars[j]]
    class.x = class(x)
    
    if(length(unique(x)) <= 2 & class.x %in% c("numeric","integer")){
      x = factor(x)
      class.x = class(x)
    }
  
    if(class.x %in% c("factor","character","logical")){
      x = factor(x)
      #gg = gam(y~x,family=fam)
      gg = gam(y~s(x,bs="re"),family=fam) 
      bp = boxplot(y~x,col='grey80',ylim=yrng,las=1)
      pp = predict(gg,data.frame(x=levels(x)),type='response') #,se.fit=TRUE)
      lines(1:nlevels(x),pp,col=2,lwd=2)  
    }
    
    if(class.x %in% c("numeric","integer")){
      uniq.x = length(unique(x))
      kmax = ifelse(uniq.x < 10, uniq.x-1, -1 ) # -1 lets gam choose kmax
      gg = gam(y~s(x,k=kmax),family=fam,xlab="")
      plot(gg,scheme=1,shift=coef(gg)[1],jit=TRUE,ylim=c(-3,3),las=1)
      abline(h=0,lty=3)
      # plot(gg,scheme=1,shift=coef(gg)[1],trans=fam$linkinv,ylim=log(range(y))) 
      # points(jitter(x),jitter(y))
    }
    title(paste(vars[j],"(",round(AIC.null-gg$aic),")"))
    AIC[j] = gg$aic
  }
  score = AIC.null - AIC  # variable importance (componentwise)
return(score)
}