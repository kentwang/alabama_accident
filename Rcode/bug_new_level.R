data <- a
fold <- cvfolds(dim(data)[1])
show.pb <- TRUE

X = model.matrix(fmla,data)
Y = as.matrix(data[,as.character(fmla[[2]])])
offset = model.offset(model.frame(fmla,data))  

fit = glm(fmla,data=a,family=poisson)

null.fmla = update(fmla,~1)  
mu = numeric(nrow(data)) 
K = sort(unique(fold))
if(show.pb) pb = txtProgressBar(style=3,min=min(K),max=max(K))
for(k in K) {
  test = which(fold == k)
  train = which(fold != k)
  
  fit0 = glm.my(fmla,data=data[train,],family=poisson)
  fit = step(fit0, trace = 0)
  mu[test] = predict(fit,newdata=data[test,],type="response")
  if(show.pb) setTxtProgressBar(pb,k)
}

Xy <- data.frame(as.matrix.data.frame(X), y=Y)
bglm <- bestglm(Xy[, -1])

#-------------------------------------------------
#- Diagnostics: Error - Factor has new levels
#  - records of '5' have been removed. And '4'?
#-------------------------------------------------
levels(a[, "LegTCType"])
which(a[, "LegTCType"] == '5') # there is only one level 5

unique(a[train, "LegTCType"]) #[1] 2 1 4
unique(a[test, "LegTCType"]) #[1] 1 2 5
