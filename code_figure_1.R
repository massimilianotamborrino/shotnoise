#----------------------------------------------------------
# Run this file to reproduce Fig. 1 of  Tamborrino and Lansky, Physica D, 418,132845, 2021 (TL2021).
# This runs a trajectory of the OU-Gamma and the OU-IG process,
# as well as of a trajectory of the shot noise process with Gamma and IG distributed jumps.

# We want to simulate Ntime points in [0,Tmax]
Tmax<-100
Ntime<- 10^5 #number of equidistant points to simulate in [0,Tmax]
x0<- 0 # non-negative starting condition of the process
alpha<-1 # positive exponential decay rate of the processes
mu<-5 #non-negative quantity entering into the jump distributions of the shot noise processes with
# Gamma and IG jumps, as well as in the Lévy measures of the OU-Gamma and OU-IG processes.
# Please refer to Table 1 of TL2021
tilde.sigma2<-3 #non-negative quantity entering into the jump distributions of the shot noise processes with
# Gamma and IG jumps, as well as in the Lévy measures of the OU-Gamma and OU-IG processes.
# Please refer to Table 1 of TL2021
scenario<- 2 # To specify that we want to simulate the shot noise process with Gamma distributed jumps

set.seed(7)
XsimG<-OU_Gamma_trajectory(Tmax,Ntime,x0,alpha,mu,tilde.sigma2,scenario)
XsimIG<-OU_IG_trajectory(Tmax,Ntime,x0,alpha,mu,tilde.sigma2)
XsimshG10<-shotnoise_JnG_trajectory(Tmax,Ntime,0,alpha,mu,tilde.sigma2,5)
XsimshIG10<-shotnoise_JnIG_trajectory(Tmax,Ntime,0,alpha,mu,tilde.sigma2,5)
dev.off()

#pdf("Figure1.pdf",width=8,height=8)
par(mfrow=c(2,2),mar=c(1.5,2,1.,1),oma=c(3,0.5,0,0))
plot(seq(0,10,by=.001),XsimG[1:10001],xlab="",ylab="",pch=".",cex.axis=2,cex.main=3)
points(seq(0,10,by=.001),XsimshG10[1:10001],pch=".",col="red")
legend("bottomright",col=c("red",1),lty=1,legend=c("shot noise","OU-Gamma"),cex=1.5)
plot(seq(0,10,by=.001),XsimIG[1:10001],xlab="time",ylab="",pch=".",cex.axis=2,cex.main=3)
legend("bottomright",col=c("red",1),lty=1,legend=c("shot noise","OU-IG"),cex=1.5)
points(seq(0,10,by=.001),XsimshIG10[1:10001],pch=".",col="red")
plot(seq(0,.8,by=.001),XsimG[1:801],xlab="",ylab="",pch=".",ylim=c(0,3.6),cex.axis=2,cex=3)
legend("bottomright",col=c("red",1),lty=1,legend=c("shot noise","OU-Gamma"),cex=1.5)
mtext("time t",side= 1,line=3,cex.axis=2,cex=2)
points(seq(0,1,by=.001),XsimshG10[1:1001],pch=".",col="red")
plot(seq(0,1,by=.001),XsimIG[1:1001],xlab="",ylab="",main="",pch=".",ylim=c(0,4),cex.axis=2,cex=3)
legend("bottomright",col=c("red",1),lty=1,legend=c("shot noise","OU-IG"),cex=1.5)
mtext("time t",side= 1,line=3,cex.axis=2,cex=2)
points(seq(0,1,by=.001),XsimshIG10[1:1001],pch=".",col="red")
#dev.off()


