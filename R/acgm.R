########################################################
# MAIN FUNCTIONS

# Part1: Estimation
# acgm.linkfunction.estim(x,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)

# Part2: Bandwidth selection under H0
# acgm.bandwidth.selection.linkfunction(data.X.H0,data.Y.H0,linkfunction.H0, kernel.function=kernel.function.epan,verbose=TRUE)

# Part3: Test statistics
# acgm.statistics(data.X,data.Y,linkfunction.H0,bandwidth='optimal',kernel.function=kernel.function.epan,integration.step=0.01 , verbose = TRUE)

# Part4: Bootstrap
# acgm.test.bootstrap(data.X,data.Y,linkfunction.H0,risk,bandwidth='optimal',kernel.function=kernel.function.epan,bootstrap=c(50,'Mammen'),integration.step=0.01,verbose=TRUE)

########################################################
# Part1: estimation

.acgm.rn0=function(x,X,h,kernel.function=kernel.function.epan)
{
  n=length(X)
  sum(kernel.function((x-X)/h))
}

.acgm.rn1=function(x,X,h,kernel.function=kernel.function.epan)
{
  n=length(X)
  sum((x-X)*kernel.function((x-X)/h))
}

.acgm.rn2=function(x,X,h,kernel.function=kernel.function.epan)
{
  n=length(X)
  sum(((x-X)**2)*kernel.function((x-X)/h))
}

.acgm.fn0=function(x,X,Y,h,kernel.function=kernel.function.epan)
{
  n=length(X)
  sum(Y*kernel.function((x-X)/h))
}

.acgm.fn1=function(x,X,Y,h,kernel.function=kernel.function.epan)
{
  n=length(X)
  sum(Y*(x-X)*kernel.function((x-X)/h))
}

acgm.linkfunction.estim=function(x,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
{
  out=rep(0,length(x))
  for (i in 1:length(x)){
    x0=x[i]
    out[i]=(.acgm.rn2(x0,data.X,bandwidth,kernel.function)*.acgm.fn0(x0,data.X,data.Y,bandwidth,kernel.function)-.acgm.rn1(x0,data.X,bandwidth,kernel.function)*.acgm.fn1(x0,data.X,data.Y,bandwidth,kernel.function))/(.acgm.rn0(x0,data.X,bandwidth,kernel.function)*.acgm.rn2(x0,data.X,bandwidth,kernel.function)-.acgm.rn1(x0,data.X,bandwidth,kernel.function)**2)
  }
  if (sum(is.na(out))>0){print('Warning: NaN in link function estimates')}
  out
}

########################################################
# Part2: bandwidth selection under H0

# Fit linkfunction under H0
.acgm.bandwidth.selection.costfunction = function( vec_bandwidth , data.X.H0 , data.Y.H0, linkfunction.H0 , kernel.function=kernel.function.epan)
{
  out=1:length(vec_bandwidth)
  compteur=0
  for (bandwidth in vec_bandwidth)
  {
    compteur=compteur+1
    x=seq( min(data.X.H0)+0.05*(max(data.X.H0)-min(data.X.H0)) , max(data.X.H0)-0.05*(max(data.X.H0)-min(data.X.H0)),length.out=100)

    err = sum((acgm.linkfunction.estim(x,data.X.H0,data.Y.H0,bandwidth,kernel.function) - linkfunction.H0(x))^2)
    if (is.na(err)){
      out[compteur] = Inf
    } else {
      out[compteur] = err
    }
  }
  out
}

acgm.bandwidth.selection.linkfunction=function(data.X.H0,data.Y.H0,linkfunction.H0, kernel.function=kernel.function.epan,verbose=TRUE)
{
  st = sort(data.X.H0,index.return=TRUE)
  dat.X.H0=st$x
  dat.Y.H0=data.Y.H0[st$ix]

  cost=function(bd)
  {
    .acgm.bandwidth.selection.costfunction(bd,dat.X.H0,dat.Y.H0,linkfunction.H0 , kernel.function)
  }
  opt = optim(mean(diff(dat.X.H0)),cost,method ="Brent",lower=min(diff(dat.X.H0))/5,upper=max(diff(dat.X.H0))*5)
  hopt = opt$par
  if (verbose==TRUE)
  {
    xgrid = seq(min(dat.X.H0),max(dat.X.H0),length.out = 100)
    ygrid = acgm.linkfunction.estim(xgrid,dat.X.H0,dat.Y.H0,hopt)
    ygrid.th = linkfunction.H0(xgrid)
    plot(xgrid,ygrid,type='l',xlab='X',ylab='Y',ylim=c(min(c(ygrid,ygrid.th)),max(c(ygrid,ygrid.th))),main='Fit with optimal bandwidth')
    lines(xgrid,ygrid.th,lty=2)
    points(dat.X.H0,dat.Y.H0,pch='+',col='gray')
  }
  hopt
}

########################################################
# Part3: test statistics

acgm.statistics=function(data.X,data.Y,linkfunction.H0,bandwidth='optimal',kernel.function=kernel.function.epan,integration.step=0.01 , verbose = TRUE)
{
  if (bandwidth=='optimal')
  {
    sd_estim = sd( data.Y-linkfunction.H0(data.X) )
    data.Y.H0 = linkfunction.H0(data.X)+rnorm(length(data.X),mean=0,sd=sd_estim)
    bdw = acgm.bandwidth.selection.linkfunction(data.X,data.Y.H0 ,linkfunction.H0 , kernel.function , verbose )
    if (verbose==TRUE)
    {
      print(paste('Optimal bandwidth under H0:',bdw,sep=' '))
    }
  } else {
    bdw=bandwidth
  }

  n=length(data.X)
  grid=seq(0,1-integration.step,by=integration.step)
  error=sum((acgm.linkfunction.estim(grid,data.X,data.Y,bdw,kernel.function)-linkfunction.H0(grid))^2)*integration.step
  n*sqrt(bdw)*error
}

########################################################
# Part3: bootstrap

.acgm.resampling=function(data.X,data.Y,bandwidth,kernel.function=kernel.function.epan,bootstrap.type='Mammen')
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
  lf=acgm.linkfunction.estim(data.X,data.X,data.Y,bandwidth,kernel.function)
  epsilon_hat=data.Y-lf
  lf+u*epsilon_hat
}


.acgm.test.stat.bootstrap=function(X,Y,newY,h,kernel.function=kernel.function.epan,integration.step)
{
  n=length(X)
  grid=seq(0,1-integration.step,by=integration.step)
  error=sum((acgm.linkfunction.estim(grid,X,Y,h,kernel.function)-acgm.linkfunction.estim(grid,X,newY,h,kernel.function))^2)*integration.step
  n*sqrt(h)*error
}

acgm.test.bootstrap=function(data.X,data.Y,linkfunction.H0,risk,bandwidth='optimal',kernel.function=kernel.function.epan,bootstrap=c(50,'Mammen'),integration.step=0.01,verbose=TRUE)
{
  if (bandwidth=='optimal')
  {
    sd_estim = sd( data.Y-linkfunction.H0(data.X) )
    data.Y.H0 = linkfunction.H0(data.X)+rnorm(length(data.X),mean=0,sd=sd_estim)
    bdw = acgm.bandwidth.selection.linkfunction(data.X,data.Y.H0 ,linkfunction.H0 , kernel.function , verbose )
    if (verbose==TRUE)
    {
      print(paste('Optimal bandwidth under H0:',bdw,sep=' '))
    }
  } else {
    bdw=bandwidth
  }

  TnB=rep(0,bootstrap[1])
  for (b in 1:bootstrap[1]){
    newY=.acgm.resampling(data.X,data.Y,bdw,kernel.function,bootstrap.type=bootstrap[2])
    TnB[b]=.acgm.test.stat.bootstrap(data.X,data.Y,newY,bdw,kernel.function,integration.step)
  }
  if (sum(is.na(TnB))>0){
    out = 'NaN in bootstrap test statistics'
    pval = 'None'
    stattest = 'None'
  } else {
    threshold=quantile(TnB,1-risk)
    stattest=acgm.statistics(data.X,data.Y,linkfunction.H0,bdw,kernel.function,integration.step,verbose=FALSE)
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
