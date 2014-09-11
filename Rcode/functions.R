## Order plots by score
## add support for other classes: Date, POSIXct, ...

## still don't like the built in gam.mgcv smooths. They are not good
#  when x data has gaps

 
library(mgcv)


component.plots <- function(fmla,data,plot=TRUE,fam = poisson(),ylim=c(-2,2),
                            max.df=6){
  
  data = model.frame(fmla,data)  # make correct transformations from fmla
  offset = as.vector(model.offset(data))
  y = model.response(data)    # response vector
  vars = attr(terms(data),"term.labels") # predictor variables
  p = length(vars)
  AIC.null = gam(y~1,family=fam,offset=offset)$aic   # AIC for baseline only model
  
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  if(plot){ 
    psqrt = sqrt(p)
    if(psqrt %% 1 == 0) xy = c(psqrt,psqrt)
    else xy = c(ceiling(psqrt),floor(psqrt))
    par(mfrow=c(xy[1],xy[2]),mar=c(2,2,1.5,1))
    ylab = "" #"log(mu/offset)-Intercept" = smooth(x)
    color = col2rgb("grey80")/255
    CI.col = rgb(color[1],color[2],color[3],alpha=.6)
    line.col = "blue"
  }  
  
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
        title(paste0(vars[j]," (",round(AIC.null-gg$aic),")"))
      }
    }
    
    if(class.x %in% c("numeric","integer")){
      nx = length(unique(x))
      #kmax = ifelse(nx < 10, nx-1, -1 ) # -1 lets gam choose kmax 
      kmax = min(nx-1,max.df-1)
      gg = gam(y~s(x,k=kmax),family=fam,offset=offset,method="REML")      
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
        title(paste0(vars[j]," (",round(AIC.null-gg$aic),")"))
      }
    }
    AIC[j] = gg$aic
  }
  score = AIC.null - AIC  # variable importance (componentwise)
  return(data.frame(variable=vars,score=round(score)))
}








