#This code contains the EM-algorithm estimation for the Hawkes process with loss covariate
rm(list=ls())
#save(seq_bus,file="E:/paper/rdata/nonstationary/seq_bus.rda")
timing_bus=cumsum(as.numeric(prc_hack_bus_arr_interval))
#timing_bus=log(timing_bus)
loss_bus=prc_hack_15_bus$Total.Records
tr_loss_bus=9.59+0.57*log(loss_bus)
####EM  alpha*t*u*exp(-beta*t)
####
source("compensator.R")
sim=gamogata(N=40,a=2,mu=30,alpha=1,beta=3)
t=sim$t
z=sim$z
#####
gam_em=function(t,z, iters, tol, mu0, alpha0, beta0){
    # Set up storage and control
  # Store iter, 3 params, log-likelihood
  N = length(t)
mat <- matrix(0, N, N)
matu<-matrix(0,N,N)
 respm<- replicate(iters, 
                     matrix(0, N,N))
 loglike=numeric(length=iters)
 mu<-numeric(length=iters)
 alpha<-numeric(length=iters)
 beta=numeric(length=iters)
 mat[,1]=1
 matu[,1]=1
 for (j in 2:N){
   mat[j,2:j]=t[j]-t[1:j-1]
   matu[j,2:j]=log((t[j]-t[1:j-1])*z[1:j-1])
 }
 matt=sum(rowSums(mat))
 matuu=sum(rowSums(matu))
for (i in 1:iters){
### E-step
  resp <- matrix(0, N, N)
  resp[, 1] = mu0
     for (j in 2:N) {    
       # z <- rgamma(N, 2, 8)
       #valid_idx = which((t[j] - t[1:j-1]) > 0) # Ensure valid values for subtraction
       resp[j, 2:j] = alpha0 * exp(-beta0 * (t[j] - t[1:j-1]))
       resp[j, 2:j] = resp[j, 2:j] * z[1:j-1] * (t[j] - t[1:j-1])
       row_sum <- rowSums(resp, na.rm = TRUE) # Avoid NA
      resp[j,] = resp [j,]/ row_sum[j]
     }
     ####### 
  respp=colSums(resp, na.rm = TRUE)[1]
   mu_new =respp / t[N]
   numlog=sum(rowSums(resp*matu, na.rm = TRUE))
    numer = sum(rowSums(resp, na.rm = TRUE))
      resmat = resp * mat
      dd = sum(rowSums(resmat, na.rm = TRUE))
     beta_new = (2*numer) / dd
    alpha_new = (4*(numer^3)) / (N*(dd^2))
#### convergence
 ###log0 ###log0
    log0=-mu0*t[N]+log(mu0)*respp-(alpha0*N)/(beta0^2)+
      log(alpha0)*numer+numlog-
      beta0*sum(rowSums(resp*mat))
    # Check if loglikelihood_old is NA
   if (is.nan(log0)) {
     stop("Initial log-likelihood computation returned NA.")
     }
    ####log_new 
    resp_new <- matrix(0, N, N)
    resp_new[, 1] = mu_new
    for (j in 2:N) {    
      # valid_idx = which((t[j] - t[1:j-1]) > 0) # Ensure valid values for subtraction
      resp_new[j, 2:j] = alpha_new * exp(-beta_new * (t[j] - t[1:j-1]))
      resp_new[j, 2:j] = resp_new[j, 2:j] * z[1:j-1] * (t[j] - t[1:j-1])
       row_sum_new <- rowSums(resp_new, na.rm = TRUE) # Avoid NA
      resp_new[j,] = resp_new[j,] / row_sum_new[j]
    }
    respp_new=colSums(resp_new, na.rm = TRUE)[1]
    numer_new = sum(rowSums(resp_new, na.rm = TRUE))
    numlog_new=sum(rowSums(resp_new*matu, na.rm = TRUE))
    ###log_new
    log_new=-mu_new*t[N]+log(mu_new)*respp_new-(alpha_new*N)/(beta_new^2)+
      log(alpha_new)*numer_new+numlog_new-
      beta_new*sum(rowSums(resp_new*mat))
    ####  # If loglik_new is NA, skip the convergence check
  if (is.nan(log_new)) {
  warning("Log-likelihood is NA. Skipping iteration.")
    next
   }
    cat(max(abs(log_new -log0)))
    if ( is.nan(max(abs(log_new -log0)))){
      warning("Log-likelihood is NA. Skipping iteration.")
      next
    } 
    if ( max(abs(log_new -log0))<tol){
      break
    } 
    #####
    respm[,,i]=resp_new
    mu[i]=mu_new
    beta[i]=beta_new
    alpha[i]=alpha_new
    loglike[i]=-log_new
    mu0=mu_new
    alpha0=alpha_new
    beta0=beta_new
}
      sim=list(mu=mu,alpha=alpha,beta=beta,loglike=loglike)
      return(sim)
 }    


#####example
gam_em(t=timing_bus,z=tr_loss_bus,
       iters=100, tol=0.001, mu0=0.0013,alpha0=0.002,beta0=0.018)

compensator_gam <- function(t,z,alpha,mu,beta) {
  n=length(t)
  M=matrix(0,n,n)
  N=matrix(0,n,n)
  y12=numeric(n)
  for (j in 2:n){
    valid_idx = which((t[j] - t[1:j-1]) > 0) # Ensure valid values for subtraction
    M[j, 1:j-1][valid_idx] = exp(-beta * (t[j] - t[1:j-1])[valid_idx])*z[1:j-1]
    N[j,1:j-1][valid_idx]=(t[j] - t[1:j-1])
    y12[j]=cumsum(z)[j]-z[j]
  } 
  T=M*N
  y11=mu*t
  y13 <- rowSums(M, na.rm = TRUE)
  y14<-rowSums(T,na.rm = TRUE)
  comp=y11+ (alpha/(beta^2))*y12-(alpha/(beta^2))*y13-(alpha/beta)*y14
  #comp_exp_z=y11+(alpha/beta)*y12-(alpha/beta)*y13
  #comp_exp=y11+(alpha*(n-1))/beta-(alpha/beta)*y13
  comp_poi=y11
  sim=list(comp=comp,comp_poi=comp_poi)
  return(sim)
}
######M
com_t=compensator_gam(t=timing_med,z=rep(1,length(timing_med)),
                      mu=0.202,alpha=1.874,beta=1.37)
####
dat=data.frame(x=seq(1,length(timing_med)-1,1),
               y=as.numeric(diff(com_t$comp)),
               y2=as.numeric(diff(com_t$comp_poi)))

####test the goodness
p <- ppoints(dat$x)
q <- quantile(dat$y,p=p)
plot(qexp(p),q, type="p",
     xlab="Theoretical Quantiles",ylab="Observed Quantiles",
     main="Q-Q Plot",
     cex.lab=1.5,cex.axis=1.5, cex.main=1.5)  
qqline(q, distribution=qexp,col="red", lty=2)