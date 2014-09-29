
library(distrom)
library(MASS)
data(fgl)
fit <- dmr(NULL, fgl[,1:9], fgl$type, lambda.min.ratio=1e-4, gamma=0, tol=1e-8, family="gaussian")

## plot the individual Poisson model fit and selection
pdf(file="fgl_coef.pdf",width=8,height=4)
par(mfrow=c(2,3),mai=c(.3,.3,.6,.2),omi=c(.4,.4,0,0))
for(j in 1:6){
	plot(fit[[j]],xlab="",ylab="",col="grey70")
	mtext(names(fit)[j],font=2,line=2,cex=.9)
	# abline(v=log(fit[[j]]$lambda[which.min(AIC(fit[[j]]))]), col="darkorange",lty=2) 
	# abline(v=log(fit[[j]]$lambda[which.min(BIC(fit[[j]]))]), col="green",lty=2) 
}
mtext("log lambda", outer=TRUE, font=3,side=1,line=1)
mtext("coefficient",outer=TRUE, font=3,side=2,line=1)
dev.off()

##### OOS experiment
lambda <- fit[[1]]$lambda
K <- 20
n <- nrow(fgl)
chunks <- round(seq.int(1,n,length=K+1))
randi <- sample.int(n)
oosdev <- rep(0,K)

getdev <- function(fit){ # multinomial deviance
			eta <- predict(fit,fgl[lo,1:9])
			ylo <- cbind(1:length(lo),fgl$type[lo])
			d <- eta[ylo] - log(rowSums(exp(eta))) 
			return(mean(-2*d)) }

for(k in 1:K){ # CV loop
	lo <- randi[chunks[k]:chunks[k+1]]
	fitk <- dmr(NULL, counts=fgl$type[-lo], fgl[-lo,1:9], gamma=0,family="gaussian",
				cores=6,lambda.start=lambda[1],lambda.min.ratio=1e-4)
	oosdev[k] <- getdev(fit)
	print(k)
}
cvm <- mean(oosdev)
cvs <- sd(oosdev)/sqrt(K-1)
cvlo <- cvm-cvs
cvhi <- cvm+cvs

library(glmnet)
fld <- c(1,rep.int(1:K,times=diff(chunks)))
fld[randi] <- fld
cvnet <- cv.glmnet(as.matrix(fgl[,1:9]), fgl$type, 
			family="multinomial",foldid=fld)

pdf(file="fgl_cv.pdf", width=7,height=4)
par(mai=c(.9,.9,.6,0))
plot(cvnet, bty="n",ylim=c(1.5,3.5),xlim=c(-10.5,2),xaxt="n")
axis(1,at=c(-10,-8,-6,-4,-2))
lines(c(.5,.5),c(cvlo,cvhi),col="grey60",lwd=2)
lines(c(0,1),c(cvlo,cvlo),col="grey60",lwd=2)
lines(c(0,1),c(cvhi,cvhi),col="grey60",lwd=2)
lines(c(.4,.6),c(cvm,cvm),col="darkblue",lwd=4)
axis(1, at=1/2,labels="DMR")
dev.off()










