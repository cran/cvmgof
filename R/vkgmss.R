########################################################
# MAIN FUNCTIONS

# Part1: Estimation
# vkgmss.linkfunction.estim(x,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
# vkgmss.sd.estim()
# vkgmss.residuals.cdf.estim()

# Part2: Bandwidth selection under H0
# vkgmss.bandwidth.selection.linkfunction(data.X.H0,data.Y.H0,linkfunction.H0, kernel.function=kernel.function.epan,verbose=TRUE)

# Part3: Test statistics
# vkgmss.statistics(data.X,data.Y,linkfunction.H0,bandwidth='optimal',kernel.function=kernel.function.epan, verbose = TRUE)

# Part4: Bootstrap
# vkgmss.test.bootstrap(data.X.H0,data.Y.H0,linkfunction.H0, kernel.function=kernel.function.epan,verbose=TRUE)

########################################################
# Part 1: Estimation

vkgmss.linkfunction.estim=function(x,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
{
  out=rep(0,length(x))
  for (k in 1:length(x)){
    x0=x[k]
    rn=sum(kernel.function((x0-data.X)/bandwidth))
    fn=sum(data.Y*kernel.function((x0-data.X)/bandwidth))
    out[k]=fn/rn
  }
  if (sum(is.na(out))>0){print('Warning: NaN in link function estimates')}
  out
}

vkgmss.sd.estim=function(x,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
{
  out=rep(0,length(x))
  for (k in 1:length(x)){
    x0=x[k]
    rn=sum(kernel.function((x0-data.X)/bandwidth))
    fn=sum(data.Y*kernel.function((x0-data.X)/bandwidth))
    fn2=sum(data.Y^2*kernel.function((x0-data.X)/bandwidth))
    out[k]=sqrt(fn2/rn-(fn/rn)^2)
  }
  if (sum(is.na(out))>0){print('Warning: NaN in standard deviation estimates')}
  out
}

.vkgmss.residuals.estim=function(data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
{
  (data.Y-vkgmss.linkfunction.estim(data.X,data.X,data.Y,bandwidth,kernel.function))/vkgmss.sd.estim(data.X,data.X,data.Y,bandwidth,kernel.function)
}

vkgmss.residuals.cdf.estim=function(u,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
{
  out=1:length(u)
  for (i in 1:length(u))
  {
    out[i]=mean(.vkgmss.residuals.estim(data.X,data.Y,bandwidth,kernel.function)<=u[i])
  }
  if (sum(is.na(out))>0){print('Warning: NaN in residual distribution estimates')}
  out
}


########################################################

# Fit linkfunction under H0
.vkgmss.bandwidth.selection.costfunction = function( vec_bandwidth , data.X.H0 , data.Y.H0, linkfunction.H0 , kernel.function=kernel.function.epan)
{
  out=1:length(vec_bandwidth)
  compteur=0
  for (bandwidth in vec_bandwidth)
  {
    compteur=compteur+1

    x=seq( min(data.X.H0)+0.05*(max(data.X.H0)-min(data.X.H0)) , max(data.X.H0)-0.05*(max(data.X.H0)-min(data.X.H0)),length.out=100)
    err = sum((vkgmss.linkfunction.estim(x,data.X.H0,data.Y.H0,bandwidth,kernel.function) - linkfunction.H0(x))^2)
    if (is.na(err)){
      out[compteur]=Inf
    } else {
      out[compteur]=err
    }
  }
  out
}

vkgmss.bandwidth.selection.linkfunction=function(data.X.H0,data.Y.H0,linkfunction.H0, kernel.function=kernel.function.epan,verbose=TRUE)
{
  st = sort(data.X.H0,index.return=TRUE)
  dat.X.H0=st$x
  dat.Y.H0=data.Y.H0[st$ix]

  cost=function(bd)
  {
    .vkgmss.bandwidth.selection.costfunction(bd,dat.X.H0,dat.Y.H0,linkfunction.H0 , kernel.function)
  }
  opt = optim(mean(diff(dat.X.H0)),cost,method ="Brent",lower=min(diff(dat.X.H0))/5,upper=max(diff(dat.X.H0))*10)
  hopt = opt$par
  if (verbose==TRUE)
  {
    xgrid = seq(min(dat.X.H0),max(dat.X.H0),length.out = 100)
    ygrid = vkgmss.linkfunction.estim(xgrid,dat.X.H0,dat.Y.H0,hopt)
    ygrid.th = linkfunction.H0(xgrid)
    plot(xgrid,ygrid,type='l',xlab='X',ylab='Y',ylim=c(min(c(ygrid,ygrid.th)),max(c(ygrid,ygrid.th))),main='Fit with optimal bandwidth')
    lines(xgrid,ygrid.th,lty=2)
    points(dat.X.H0,dat.Y.H0,pch='+',col='gray')
  }
  hopt
}

########################################################

.vkgmss.residuals.H0.estim=function(data.X,data.Y,linkfunction.H0,bandwidth,kernel.function=kernel.function.epan)
{
  (data.Y-linkfunction.H0(data.X))/vkgmss.sd.estim(data.X,data.X,data.Y,bandwidth,kernel.function)
}

.vkgmss.residuals.H0.cdf.estim=function(u,data.X,data.Y,linkfunction.H0,bandwidth,kernel.function=kernel.function.epan)
{
  out=1:length(u)
  for (i in 1:length(u))
  {
    out[i]=mean(.vkgmss.residuals.H0.estim(data.X,data.Y,linkfunction.H0,bandwidth,kernel.function)<=u[i])
  }
  out
}

vkgmss.statistics=function(data.X,data.Y,linkfunction.H0,bandwidth='optimal',kernel.function=kernel.function.epan,verbose=TRUE)
{
  if (bandwidth=='optimal')
  {
    sd_estim = sd( data.Y-linkfunction.H0(data.X) )
    data.Y.H0 = linkfunction.H0(data.X)+rnorm(length(data.X),mean=0,sd=sd_estim)
    bdw = vkgmss.bandwidth.selection.linkfunction(data.X,data.Y.H0 ,linkfunction.H0 , kernel.function , verbose )
    if (verbose==TRUE)
    {
      print(paste('Optimal bandwidth under H0:',bdw,sep=' '))
    }
  } else {
    bdw=bandwidth
  }
  res.H0 = .vkgmss.residuals.H0.estim(data.X ,data.Y,linkfunction.H0,bdw,kernel.function)
  error = (vkgmss.residuals.cdf.estim(res.H0,data.X,data.Y,bdw,kernel.function) - .vkgmss.residuals.H0.cdf.estim(res.H0,data.X,data.Y,linkfunction.H0,bdw,kernel.function))^2
  sum(error)
}

########################################################

.vkgmss.resampling=function(data.X,data.Y,bandwidth,kernel.function=kernel.function.epan,bootstrap.type='Mammen')
{
  n=length(data.X)
  if (bootstrap.type=='Gaussian'){u=rnorm(n,mean=0,sd=1)}
  else if (bootstrap.type=='Mammen'){
    u=rep(0,n)
    for (k in 1:n){
      v=runif(1,min=0,max=1)
      if (v<(sqrt(5)+1)/(2*sqrt(5))){uu=-(sqrt(5)-1)/2} else {uu=(sqrt(5)+1)/2}
      u[k]=uu
    }
  }
  else if (bootstrap.type=='Rademacher'){
    u=rep(0,n)
    for (k in 1:n){
      v=runif(1,min=0,max=1)
      if (v<0.5){uu=-1} else {uu=1}
      u[k]=uu
    }
  }
  lf=vkgmss.linkfunction.estim(data.X,data.X,data.Y,bandwidth,kernel.function)
  epsilon_hat=data.Y-lf
  lf+u*epsilon_hat
}


.vkgmss.test.stat.bootstrap=function(data.X,data.Y,data.newY,bandwidth,kernel.function=kernel.function.epan,integration.step)
{
  res.H0 = .vkgmss.residuals.estim(data.X,data.newY,bandwidth,kernel.function)
  error = (vkgmss.residuals.cdf.estim(res.H0,data.X,data.Y,bandwidth,kernel.function) - vkgmss.residuals.cdf.estim(res.H0,data.X,data.newY,bandwidth,kernel.function))^2
  sum(error)
}

vkgmss.test.bootstrap=function(data.X,data.Y,linkfunction.H0,risk,bandwidth='optimal',kernel.function=kernel.function.epan,bootstrap=c(50,'Mammen'),verbose=TRUE)
{
  if (bandwidth=='optimal')
  {
    sd_estim = sd( data.Y-linkfunction.H0(data.X) )
    data.Y.H0 = linkfunction.H0(data.X)+rnorm(length(data.X),mean=0,sd=sd_estim)
    bdw = vkgmss.bandwidth.selection.linkfunction(data.X,data.Y.H0 ,linkfunction.H0 , kernel.function , verbose )
    if (verbose==TRUE)
    {
      print(paste('Optimal bandwidth under H0:',bdw,sep=' '))
    }
  } else {
    bdw=bandwidth
  }

  TnB=rep(0,bootstrap[1])
  for (b in 1:bootstrap[1]){
    newY=.vkgmss.resampling(data.X,data.Y,bdw,kernel.function,bootstrap.type=bootstrap[2])
    TnB[b]=.vkgmss.test.stat.bootstrap(data.X,data.Y,newY,bdw,kernel.function)
  }
  if (sum(is.na(TnB))>0){
    out = 'NaN in bootstrap test statistics'
    pval = 'None'
    stattest = 'None'
  } else {
    threshold=quantile(TnB,1-risk)
    stattest=vkgmss.statistics(data.X,data.Y,linkfunction.H0,bdw,kernel.function,verbose=FALSE)
    if (is.na(stattest)==TRUE){
      out='test statistics is NaN'
      pval='None'
    } else {
      if(stattest<threshold){out='accept H0'} else {out='reject H0'}
      pval=sum(TnB>stattest)/length(TnB)
      if (pval==0){
        pval=paste('<',1/length(TnB),sep=' ')
      }
    }
  }
  list('decision'=out,'bandwidth'=bdw,'pvalue'=pval,'test_statistics'=stattest)
}
########################################################
