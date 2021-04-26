#----------------------------------------------------------
# Run this file to reproduce Fig. 2 of  Tamborrino and Lansky, Physica D, 418,132845, 2021 (TL2021).
# This plots the probability density functions of the shot noise processes Xn, the Gaussian OU processes and the
# the non-Gaussian OU processes (10 of TL2021) at time t=T0=1 and t=T0b=100, when the underlying jumps of the shot noise process are
# Bernoulli/Poissonian, Gamma or Inverse-Gaussian distributed. The non-Gaussian OU-processes are then
# OU-Poisson, OU-Gamma and OU-IG, respectively

# We want to simulate Ntime points in [0,Tmax]
Tmax<-100
Ntime<- 10^5 #number of equidistant points to simulate in [0,Tmax]
alpha<-1
mu<-5
scenario<- 2 # To specify that we want to simulate the shot noise process with Gamma distributed jumps

lambda0<-1000 #positive intensity of the underlying Poisson process of the shot noise processes
alpha0<-1# positive exponential decay rate of the processes
mu0<-2 #non-negative quantity entering into the jump distributions of the shot noise processes with, as well as in the Gaussian OU and
# in the Lévy measures of the non-Gaussian OU processes
x0<- 0 # non-negative starting condition of the process
tilde.sigma2<-3 #non-negative quantity entering into the jump distributions of the shot noise processes with, as well as in the Gaussian OU and
# in the Lévy measures of the non-Gaussian OU processes
# Please refer to Table 1 of TL2021
T0<-1 #Time at which we want to evaluate the densities of the processes.
T0b<-100 #Other time at which we want to evaluate the densities of the processes, with the idea of being in the asymptotic regimes.
Nsim<-10^6 # Number of simulation for each process
Npoints<-10^3 #number of equidistant points to simulate in [0,t] for the OU-Gamma and OU-IG processes

# shotnoise_JnBer code to simulate from a shot noise process with Bernoulli distributed jumps
# shotnoise_JnP code to simulate from a shot noise process with Poissonian distributed jumps
# shotnoise_JnG code to simulate from a shot noise process with Gamma distributed jumps
# shotnoise_JnIG code to simulate from a shot noise process with inverse-Gaussian distributed jumps
# OU_Poisson code to simulate from a non-Gaussian OU-Poisson process
# OU_Gamma code to simulate from a non-Gaussian OU-Gamma process
# OU_IG code to simulate from a non-Gaussian OU-IG process

set.seed(7)
D_BEx00 <-Vectorize(shotnoise_JnBer)(T0,Nsim,x0,alpha0,mu0,lambda0)
D_Px00 <-Vectorize(shotnoise_JnP)(T0,Nsim,x0,alpha0,mu0,lambda0)
D_Gx00 <-Vectorize(shotnoise_JnG)(T0,Nsim,x0,alpha0,mu0,tilde.sigma20,lambda0)
D_IGx00 <- Vectorize(shotnoise_JnIG)(T0,Nsim,x0,alpha0,mu0,tilde.sigma20,lambda0)

D_BEsimx00 <-Vectorize(OU_Poisson)(T0,Nsim,x0,alpha0,mu0)
XOUGx00 <-Vectorize(OU_Gamma)(T0,Npoints,Nsim,x0,alpha0,mu0,tilde.sigma20,2)
XOUIGx00 <-Vectorize(OU_IG)(T0,Npoints,Nsim,x0,alpha0,mu0,tilde.sigma20)
D_BEb <-Vectorize(shotnoise_JnBer)(T0b,Nsim,x0,alpha0,mu0,lambda0)
D_Pb<-Vectorize(shotnoise_JnP)(T0b,Nsim,x0,alpha0,mu0,lambda0)

D_Gb <-Vectorize(shotnoise_JnG)(T0b,Nsim,x0,alpha0,mu0,tilde.sigma20,lambda0)
D_IGb <-Vectorize(shotnoise_JnIG)(T0b,Nsim,x0,alpha0,mu0,tilde.sigma20,lambda0)

D_BEsimb <-Vectorize(OU_Poisson)(T0b,Nsim,x0,alpha0,mu0)
XOUGb <-Vectorize( OU_Gamma)(T0b,Npoints,Nsim,x0,alpha0,mu0,tilde.sigma20,2)
XOUIGb <-Vectorize(OU_IG)(T0b,Npoints,Nsim,x0,alpha0,mu0,tilde.sigma20)

#pdf("Figure2.pdf",width=12,height=8)
par(mfrow=c(2,3),mar=c(2,3,2.5,1),oma=c(0,2.5,0,0))
plot(density(D_BEsimx00),col=1,type="l",xlim=c(0,6),ylim=c(0,1),cex=3,xlab="",ylab="",
     main=expression(J[n]:~Bernoulli / Poisson),lwd=3,cex.axis=3,cex.main=3.3)
mtext(text="t=1",outer=F,side=2,line=3,cex=2)
lines(density(D_BEx00),col="grey",lwd=1.5)
curve(dnorm(x,x0*exp(-alpha0*T0)+mu0*(1-exp(-alpha0*T0))/alpha0, sd=sqrt(mu0*(1-exp(-2*alpha0*T0))/(2*alpha0))),add=T,col="blue")
legend("topright",lty=1,col=c(1,"grey","blue"),legend=c("shot noise","OU-Poisson","OU"),cex=2.1)

plot(density(D_Gx00),col=1,lwd=3,xlim=c(0,6),ylim=c(0,.95),type="l",cex=3,xlab="",ylab="",
     main=expression(J[n]:~ Gamma),cex.axis=3,cex.main=3.3)
lines(density(XOUGx00),col="grey",lwd=1.5)
curve(dnorm(x,x0*exp(-alpha0*T0)+mu0*(1-exp(-alpha0*T0))/alpha0, sd=sqrt(tilde.sigma20*(1-exp(-2*alpha0*T0))/(2*alpha0))),add=T,col="blue")
legend("topright",lty=1,col=c(1,"grey","blue"),legend=c("shot noise","OU-Gamma","OU"),cex=2.1)

plot(density(D_IGx00),col=1,type="l",xlim=c(0,6),ylim=c(0,.95),cex=3,xlab="",ylab="",
     main=expression(J[n]:~ IG),cex.axis=3,lwd=3,cex.main=3.3)
lines(density(XOUIGx00),col="grey")
curve(dnorm(x,x0*exp(-alpha0*T0)+mu0*(1-exp(-alpha0*T0))/alpha0, sd=sqrt(tilde.sigma20*(1-exp(-2*alpha0*T0))/(2*alpha0))),add=T,col="blue")
legend("topright",lty=1,col=c(1,"grey","blue"),legend=c("shot noise","OU-IG","OU"),cex=2.1)

plot(density(D_BEb),col=1,lty=1,type="l",xlim=c(0,6),ylim=c(0,.55),
     main="",xlab="",ylab="",cex=3,cex.axis=3,lwd=3)
mtext(text="t=100",outer=F,side=2,line=3,cex=2)
lines(density(D_BEsimb),col="grey",lty=1)
curve(dnorm(x,x0*exp(-alpha0*T0)+mu0*(1-exp(-alpha0*T0b))/alpha0, sd=sqrt(mu0*(1-exp(-2*alpha0*T0b))/(2*alpha0))),add=T,col="blue")
legend("topright",lty=1,col=c(1,"grey","blue"),legend=c("shot noise","OU-Poisson","OU"),cex=2.)

plot(density(D_Gb),col=1,lwd=3,xlim=c(0,6),ylim=c(0,.55),type="l",cex=3,xlab="",ylab="",
     main="",cex.axis=3)
lines(density(XOUGb),col="grey")
curve(dnorm(x,x0*exp(-alpha0*T0)+mu0*(1-exp(-alpha0*T0b))/alpha0, sd=sqrt(tilde.sigma20*(1-exp(-2*alpha0*T0b))/(2*alpha0))),add=T,col="blue")
legend("topright",lty=1,col=c(1,"grey","blue"),legend=c("shot noise","OU-Gamma","OU"),cex=2.1)

plot(density(D_IGb),col=1,type="l",xlim=c(0,6),ylim=c(0,.55),cex=3,xlab="",ylab="",
     main="",cex.axis=3,lwd=3,)
lines(density(XOUIGb),col="grey")
curve(dnorm(x,x0*exp(-alpha0*T0)+mu0*(1-exp(-alpha0*T0b))/alpha0, sd=sqrt(tilde.sigma20*(1-exp(-2*alpha0*T0b))/(2*alpha0))),add=T,col="blue")
legend("topright",lty=1,col=c(1,"grey","blue"),legend=c("shot noise","OU-IG","OU"),cex=2.1)
#dev.off()

