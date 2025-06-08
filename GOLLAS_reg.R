#  beta = ( -log(2)/log(2/pi*asin(sqrt(mu))) )
# alpha = sigma
lgollas<-expression( 
  log(
  ( sigma * ( -log(2)/log(2/pi*asin(sqrt(mu))) ) * ( (2 / pi) * asin(sqrt(y)) )^(sigma * ( -log(2)/log(2/pi*asin(sqrt(mu))) ) - 1) * 
      (1 - ( (2 / pi) * asin(sqrt(y)) )^( -log(2)/log(2/pi*asin(sqrt(mu))) ))^(sigma - 1)  )/
       (  pi * sqrt( y - y^2 ) * ( ( (2 / pi) * asin(sqrt(y)) )^(sigma * ( -log(2)/log(2/pi*asin(sqrt(mu))) )) + 
                                              (1 - ( (2 / pi) * asin(sqrt(y)) )^( -log(2)/log(2/pi*asin(sqrt(mu))) ))^sigma )^2 )
    ) # close log 
) # close expression


m1<-D(lgollas,"mu")
s1<-D(lgollas,"sigma")
ms2<-D(m1,"sigma")
m1

GOLLAS<-function (mu.link = "logit", sigma.link = "log") 
{
  mstats <- checklink("mu.link", "GOLLAS", substitute(mu.link), 
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "GOLLAS", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("GOLLAS", "lgollas"), 
                 parameters = list(mu = TRUE, sigma = TRUE), 
                 nopar = 2, 
                 type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 
                 dldm = function(y, mu, sigma) {#ok
                   a<-0
                   b<-1
                   dldm <- eval(m1)
                   dldm
                 }, 
                 d2ldm2 = function(y,mu, sigma) {
                   a<-0
                   b<-1
                   dldm <- eval(m1)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                   d2ldm2
                 }, 
                 dldd = function(y, mu, sigma) {#ok
                   a<-0
                   b<-1
                   dldd <- eval(s1)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   a<-0
                   b<-1
                   dldd <- eval(s1)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   a<-0
                   b<-1
                   dldm <- eval(m1)
                   dldd <- eval(s1)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd  
                 }, 
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dGOLLAS(y=y, mu=mu, sigma=sigma)), 
                 rqres = expression(
                   rqres(pfun = "pGOLLAS",  type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 
                 mu.initial = expression(     mu <- rep(median(y),length(y))),   
                 sigma.initial = expression(sigma<- rep(.5, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y > 0 &  y < 1)
  ), 
  class = c("gamlss.family", "family"))
}


# density function
dGOLLAS<-function(y, mu = 0.7, sigma = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  fy1 <- ( sigma * ( -log(2)/log(2/pi*asin(sqrt(mu))) ) * ( (2 / pi) * asin(sqrt(y)) )^(sigma * ( -log(2)/log(2/pi*asin(sqrt(mu))) ) - 1) * 
             (1 - ( (2 / pi) * asin(sqrt(y)) )^( -log(2)/log(2/pi*asin(sqrt(mu))) ))^(sigma - 1)  )/
    (  pi * sqrt( y - y^2 ) * ( ( (2 / pi) * asin(sqrt(y)) )^(sigma * ( -log(2)/log(2/pi*asin(sqrt(mu))) )) + 
                                  (1 - ( (2 / pi) * asin(sqrt(y)) )^( -log(2)/log(2/pi*asin(sqrt(mu))) ))^sigma )^2 )
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}
#------------------------------------------------------------------------------------------ #ok
# cumulative distribution function
pGOLLAS<-function(q, mu = 0.7, sigma = 0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))

  cdf1 <- ( ( (2 / pi) * asin(sqrt(y)) )^(sigma * ( -log(2)/log(2/pi*asin(sqrt(mu))) )) )/
    ( ( (2 / pi) * asin(sqrt(y)) )^(sigma * ( -log(2)/log(2/pi*asin(sqrt(mu))) ))+
        (1-( (2 / pi) * asin(sqrt(y)) )^( -log(2)/log(2/pi*asin(sqrt(mu))) ))^sigma )
  
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#------------------------------------------------------------------------------------------ #ok
# quantile function
qGOLLAS<-function(u,mu,sigma)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(u <= 0) | any(u >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))

    q <- sin(pi * (  (1 / (( ( (1 - u) / u )^(1 / sigma) ) + 1))^(1 / ( -log(2)/log(2/pi*asin(sqrt(mu))) ))  ) / 2)^2
    q
}

# inversion method for randon generation
rGOLLAS<-function(n,mu,sigma)
{
  u<- runif(n)
  y<- sin(pi * (  (1 / (( ( (1 - u) / u )^(1 / sigma) ) + 1))^(1 / ( -log(2)/log(2/pi*asin(sqrt(mu))) ))  ) / 2)^2
  y
}

# set.seed(100)
# library(gamlss)
# n<-1000
# X<-runif(n)
# ddd<-make.link("logit")
# y<-rGOLLAS(n,ddd$linkinv(-1.6+1.25*X),2.3)
# 
# saida<-gamlss(y~X, family="GOLLAS")
# sai<-summary(saida)
# round(sai,4)
