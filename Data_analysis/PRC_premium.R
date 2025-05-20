# This code is used to compute the premiums under two different principles
rm(list=ls())
library(hawkesbow)
library(ggplot2)
#max_index <- which.max(prc_hack_15_med$Total.Records)
#prc_hack_15_med <-prc_hack_15_med[-max_index,]
timing_bus=cumsum(as.numeric(prc_hack_bus_arr_interval))
loss_exp=function(mu,mu_z,alpha,beta,t){
  kappa=beta-alpha*mu_z*mu
  loss=(mu*mu_z+alpha*(mu_z^2)/kappa)*t-alpha*(mu_z^2)/(kappa^2)
  return(loss)
}
#### alpha,mu,beta are different from the exponential
loss_gam=function(mu,mu_z,alpha,beta,t){
  kappa=beta-alpha*mu_z*mu
  los=alpha*mu*(mu_z^2)/(beta*kappa)
  lossum=(1/kappa)+(1/beta)+(1/(beta+kappa))
  loss=(mu*mu_z+los)*t-los*lossum
  return(loss)
}
######
haw<-hawkes_ogata(1440,0.28,0.0025,17)
dat<-haw$p
#####lossexp
losexpz1=loss_exp(0.042,11.34,0.014,4.67,t=timing_bus)
losexpz2=loss_exp(0.077,18.23,0.018,7.5,t=timing_bus)
losexpz3=loss_exp(0.159,24.27,0.026,10.8,t=timing_bus)
x1=rep(timing_bus,times=3)
group=rep(c("q=0.25","q=0.5","q=0.75"),each=length(timing_bus))
y=c(losexpz1,losexpz2,losexpz3)
lossexp=data.frame(x1,group,y)
lossexp2=data.frame(rep(timing_bus,times=2),
                    group=rep(c("Exponential kernel","General kernel"),
                              each=length(timing_bus)),
                    losexpz1,losexpz2,losexpz3)
#####
f=ggplot(lossexp,aes(x=x1,y=y,linetype=group))+geom_step()+geom_point()+
  theme_classic()
f=f+theme(axis.title = element_text(size = 15))+
  theme(axis.text = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))
f=f+scale_linetype_manual(values=c("solid","dotted","dashed"))
#f=f+scale_shape_discrete(labels=c("α=0.2","α=0.3","α=0.4","α=0.5"))
f=f+theme(legend.title = element_blank())
f=f+labs(y="Aggregate Loss",x="Time")
f=f+theme(legend.position = "bottom")
#######lossgam
losgamz1=loss_gam(0.20,12.34,0.024,8.67,t=timing_bus)
losgamz2=loss_gam(0.31,15.23,0.062,16.5,t=timing_bus)
losgamz3=loss_gam(0.59,18.27,0.082,20.8,t=timing_bus)
x11=rep(timing_bus,times=3)
group1=rep(c("lambda1","lambda2","lambda3"),each=length(timing_bus))
y1=c(losgamz1,losgamz2,losgamz3)
lossgam=data.frame(x11,group1,y1)
lossgam2=data.frame(seq(1:length(timing_bus)),losgamz1,losgamz2,losgamz3)
#####plotlossgam
###### data 
re=rep(c("lambda1","lambda2","lambda3"),
       each=length(timing_bus))
lossexp2=data.frame(rep(seq(1,length(timing_bus)),times=6),
  group=rep(c("Exponential kernel","General kernel"),each=3
            *length(timing_bus)),
      y= c(losexpz1,losexpz2,losexpz3,losgamz1,losgamz2,losgamz3),
repp=rep(re,times=2))
#####box plot
ggplot(lossexp2, aes(x = repp, y = y, fill = group)) + 
  geom_boxplot()+ theme_classic()+
  theme(legend.title = element_blank())+
  theme(axis.title = element_text(size = 15))+
  theme(axis.text = element_text(size = 15))+
  theme(legend.text = element_text(size = 15))+
  scale_fill_brewer(palette = "Pastel2") + #利用RColorBrewer中Pastel2配色对分组变量进行填充
  labs(x= "Quantile points", y = "Premiums") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") #legend.position用来设置图例放置的位置（绘图区域内）        