# An comparison of two kinds of kernel function
library(latex2exp)
fun1<-function(x){
  alpha=0.5
  beta=0.8
  fun1=alpha*exp(-beta*x)
}
  fun2<-function(x){
    alpha=0.5
    beta=0.8
   fun2=alpha*x*exp(-beta*x)
  }
  curve(expr=fun1,0,8,xlab="Time",ylab="Decay kernel",cex.lab=2,
        cex.axis=2,lwd=1,col="red",lty=2,mgp=c(2,1,0))
  curve(expr=fun2,0,8,add=TRUE,col="black",lty=1)
  #abline(h=1,col="blue",lty=3)
  legend('topright', 
         legend=TeX(c('$\\phi(a)=\\alpha a e^{-\\beta a}$','$\\phi(a)=\\alpha e^{-\\beta a}$')), lty=c(1, 2), 
         col=c("black","red"), cex=2,bty="n",
         ncol=2, inset=0.001)