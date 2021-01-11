# Kernel functions

kernel.function.epan=function(u)
{
  3/4*(1-u^2)*(abs(u)<=1)
}

kernel.function.gauss=function(u)
{
  dnorm(u)
}

kernel.function.quart=function(u)
{
  15/16*(1-u^2)*(abs(u)<=1)
}
