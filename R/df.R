########################################################
# MAIN FUNCTIONS

# Part1: Estimation
# df.cdf.estim(x,y,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
# df.linkfunction.estim(x,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)

# Part2: Bandwidth selection under H0
# df.bandwidth.selection.linkfunction(data.X.H0,data.Y.H0,linkfunction.H0 , kernel.function=kernel.function.epan,integration.step=0.01,verbose=TRUE)

# Part3: Test statistics
# df.statistics(data.X,data.Y,cdf.H0,bandwidth,kernel.function=kernel.function.epan,integration.step=0.01)

# Part4: Bootstrap
# df.test.bootstrap(data.X,data.Y,cdf.H0,risk,bandwidth,kernel.function=kernel.function.epan,bootstrap=c(50,'Mammen'),integration.step=0.01)

########################################################
# Part1: estimation

.df.rn0=function(x,X,h,kernel.function=kernel.function.epan)
{
  mat=outer(x,X,'-')
  out=apply(kernel.function(mat/h),1,sum,na.rm=TRUE)
  out
}

.df.rn1=function(x,X,h,kernel.function=kernel.function.epan)
{
  mat=outer(x,X,'-')
  out=apply(mat*kernel.function(mat/h), 1,sum, na.rm=TRUE)
  out
}

.df.rn2=function(x,X,h,kernel.function=kernel.function.epan)
{
  mat=outer(x,X,'-')
  out=apply(mat^2*kernel.function(mat/h), 1,sum, na.rm=TRUE)
  out
}

.df.fn0=function(x,y,X,Y,h,kernel.function=kernel.function.epan)
{
  mat=outer(x,X,'-')
  out=matrix(0,nrow=length(x),ncol=length(y))
  for (k in 1:length(y)){out[,k]=kernel.function(mat/h)%*%(Y<=y[k])}
  out
}

.df.fn1=function(x,y,X,Y,h,kernel.function=kernel.function.epan)
{
  mat=outer(x,X,'-')
  out=matrix(0,nrow=length(x),ncol=length(y))
  for (k in 1:length(y)){out[,k]=(mat*kernel.function(mat/h))%*%(Y<=y[k])}
  out
}

df.cdf.estim=function(x,y,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
{
  out=(.df.rn2(x,data.X,bandwidth)*.df.fn0(x,y,data.X,data.Y,bandwidth)-.df.rn1(x,data.X,bandwidth)*.df.fn1(x,y,data.X,data.Y,bandwidth))/(.df.rn2(x,data.X,bandwidth)*.df.rn0(x,data.X,bandwidth)-.df.rn1(x,data.X,bandwidth)^2)
  if (sum(is.na(out))>0){print('Warning: NaN in cumulative distribution function estimates')}
  out
}

df.linkfunction.estim=function(x,data.X,data.Y,bandwidth,kernel.function=kernel.function.epan)
{
  n=length(data.X)
  sortY=sort(data.Y)
  Fdrestim=df.cdf.estim(x,sortY,data.X,data.Y,bandwidth,kernel.function)
  diffsortY=sortY[2:n]-sortY[1:(n-1)]
  out=sortY[n]*Fdrestim[,n] - sortY[1]*Fdrestim[,1] - Fdrestim[,1:(n-1)]%*%diffsortY
  if (sum(is.na(out))>0){print('Warning: NaN in link function estimates')}
  out
}


########################################################
# Part2: bandwidth selection under H0

# Fit linkfunction under H0
.df.bandwidth.selection.costfunction.linkfunction = function( vec_bandwidth , data.X.H0,data.Y.H0,linkfunction.H0 , kernel.function=kernel.function.epan)
{
  out=1:length(vec_bandwidth)
  compteur=0
  for (bandwidth in vec_bandwidth)
  {
    compteur=compteur+1
    x=seq( min(data.X.H0)+0.05*(max(data.X.H0)-min(data.X.H0)) , max(data.X.H0)-0.05*(max(data.X.H0)-min(data.X.H0)),length.out=100)
    err = sum((df.linkfunction.estim(x,data.X.H0,data.Y.H0,bandwidth,kernel.function) - linkfunction.H0(x))^2)
    if (is.na(err)){
      out[compteur]=Inf
    } else {
      out[compteur]=err
    }
  }
  out
}

df.bandwidth.selection.linkfunction=function(data.X.H0,data.Y.H0,linkfunction.H0, kernel.function=kernel.function.epan,verbose=TRUE)
{
  st = sort(data.X.H0,index.return=TRUE)
  dat.X.H0=st$x
  dat.Y.H0=data.Y.H0[st$ix]

  cost=function(bd)
  {
    .df.bandwidth.selection.costfunction.linkfunction(bd,dat.X.H0,dat.Y.H0,linkfunction.H0 , kernel.function)
  }
  opt = optim(mean(diff(dat.X.H0)),cost,method ="Brent",lower=min(diff(dat.X.H0))/5,upper=max(diff(dat.X.H0))*5)
  hopt = opt$par
  if (verbose==TRUE)
  {
    xgrid = seq(min(dat.X.H0),max(dat.X.H0),length.out = 100)
    ygrid = df.linkfunction.estim(xgrid,dat.X.H0,dat.Y.H0,hopt)
    ygrid.th = linkfunction.H0(xgrid)
    plot(xgrid,ygrid,type='l',xlab='X',ylab='Y',ylim=c(min(c(ygrid,ygrid.th)),max(c(ygrid,ygrid.th))),main='Fit with optimal bandwidth')
    lines(xgrid,ygrid.th,lty=2)
    points(dat.X.H0,dat.Y.H0,pch='+',col='gray')
  }
  hopt
}


########################################################
# Part3: test statistics

.df.tnx=function(x,X,Y,h,cdf.H0,kernel.function=kernel.function.epan)
{
  n=length(X)
  sortY=sort(Y)
  Fdr0mat=cdf.H0(x,sortY)

  #First term of Tnx
  Fdrestim2=df.cdf.estim(x,sortY,X,Y,h,kernel.function)^2
  etape1=Fdr0mat[,2:n]-Fdr0mat[,1:(n-1)]
  etape2=Fdrestim2[,1:(n-1)]*etape1
  etape3=apply(etape2,1,sum)

  #Second term of Tnx
  Fdr0mat2=Fdr0mat^2
  etape4=Fdr0mat2[,2:n]-Fdr0mat2[,1:(n-1)]
  etape5=df.cdf.estim(x,sortY,X,Y,h,kernel.function)[,1:(n-1)]*etape4
  etape6=apply(etape5,1,sum)

  #Third term of Tnx
  etape7=Fdr0mat2[,n]-Fdr0mat[,n]

  1/3+etape3-etape6+etape7
}

df.statistics=function(data.X,data.Y,cdf.H0,bandwidth,kernel.function=kernel.function.epan,integration.step=0.01)
{
  n=length(data.X)
  grid=seq(0,1-integration.step,by=integration.step)
  error=sum(.df.tnx(grid,data.X,data.Y,bandwidth,cdf.H0,kernel.function))*integration.step
  n*sqrt(bandwidth)*error
}

########################################################
# Part4: bootstrap

.df.resampling=function(data.X,data.Y,bandwidth,kernel.function=kernel.function.epan,bootstrap.type='Mammen')
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
  lf=df.linkfunction.estim(data.X,data.X,data.Y,bandwidth,kernel.function)
  epsilon_hat=data.Y-lf
  lf+u*epsilon_hat
}

.df.integration_Fn=function(matrice_a_integrer,X,x,h,kernel.function=kernel.function.epan)
{
  n=length(X)
  mat=outer(x,X,'-')
  matK=kernel.function(mat/h)
  apply(matrice_a_integrer*(.df.rn2(x,X,h,kernel.function)*matK-.df.rn1(x,X,h,kernel.function)*matK*mat)/(.df.rn2(x,X,h,kernel.function)*.df.rn0(x,X,h,kernel.function)-.df.rn1(x,X,h,kernel.function)^2),1,sum)
}

.df.tnxb=function(x,newY,X,Y,h,kernel.function=kernel.function.epan,integration.step=0.01)
{
  n=length(X)
  Fdrestim=df.cdf.estim(x,newY,X,Y,h,kernel.function)
  Fdrestimboot=df.cdf.estim(x,newY,X,newY,h,kernel.function)

  mat_a_int=(Fdrestim-Fdrestimboot)^2
  .df.integration_Fn(mat_a_int,X,x,h)
}

.df.test.stat.bootstrap=function(X,Y,newY,h,kernel.function=kernel.function.epan,integration.step=0.01)
{
  n=length(X)
  grid=seq(0,1-integration.step,by=integration.step)
  error=sum(.df.tnxb(grid,newY,X,Y,h,kernel.function,integration.step))*integration.step
  n*sqrt(h)*error
}

df.test.bootstrap=function(data.X,data.Y,cdf.H0,risk,bandwidth,kernel.function=kernel.function.epan,bootstrap=c(50,'Mammen'),integration.step=0.01)
{
  TnB=rep(0,bootstrap[1])
  for (b in 1:bootstrap[1]){
    newY=.df.resampling(data.X,data.Y,bandwidth,kernel.function,bootstrap.type=bootstrap[2])
    TnB[b]=.df.test.stat.bootstrap(data.X,data.Y,newY,bandwidth,kernel.function,integration.step)
  }
  if (sum(is.na(TnB))>0){
    out = 'NaN in bootstrap test statistics'
    pval = 'None'
    stattest = 'None'
  } else {
    threshold=quantile(TnB,1-risk)
    stattest=df.statistics(data.X,data.Y,cdf.H0,bandwidth,kernel.function,integration.step)
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
  list('decision'=out,'bandwidth'=bandwidth,'pvalue'=pval,'test_statistics'=stattest)
}
########################################################
