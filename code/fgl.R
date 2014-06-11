
library(dmr)
library(MASS)
data(fgl)
fit <- dmr(NULL, fgl[,1:9], fgl$type, store=TRUE, lambda.min.ratio=1e-4, gamma=0)
## do grouped BIC model selection
B <- coef(fit, grouped=TRUE)
log(B@lambda)
Bbic <- coef(fit, grouped=TRUE, k=log(nrow(fgl)))

## plot the individual Poisson model fit and selection
pdf(file="fgl_coef.pdf",width=8,height=4)
par(mfrow=c(2,3),mai=c(.3,.3,.6,.2),omi=c(.4,.4,0,0))
for(j in 1:6){
	plot(fit[[j]],xlab="",ylab="",select=FALSE)
	mtext(names(fit)[j],font=2,line=2,cex=.9)
	abline(v=log(fit[[j]]$lambda[which.min(AIC(fit[[j]]))]), col="darkorange",lty=2) 
	abline(v=log(fit[[j]]$lambda[which.min(BIC(fit[[j]]))]), col="green",lty=2) 
}
mtext("log lambda", outer=TRUE, font=3,side=1,line=1)
mtext("coefficient",outer=TRUE, font=3,side=2,line=1)
dev.off()

lambda <- fit[[1]]$lambda
K <- 20
n <- nrow(fgl)
chunks <- round(seq.int(1,n,length=K+1))
D <- matrix(nrow=K, ncol=100)
randi <- sample.int(n)
icD <- matrix(nrow=K,ncol=2)
colnames(icD) <- c("AIC","BIC")

getdev <- function(s){
			eta <- predict(fitk,fgl[lo,1:9],select=s)
			ylo <- cbind(1:length(lo),fgl$type[lo])
			d <- eta[ylo] - log(rowSums(exp(eta))) 
			return(mean(-2*d)) }

for(k in 1:K){
	lo <- randi[chunks[k]:chunks[k+1]]
	fitk <- dmr(counts=fgl$type[-lo], fgl[-lo,1:9], gamma=0,
				cores=6,lambda.start=lambda[1],lambda.min.ratio=1e-4)
	D[k,] <- sapply(1:100, getdev)
	icD[k,'AIC'] <- getdev(sapply(fitk,function(m) which.min(AIC(m))))
    icD[k,'BIC'] <- getdev(sapply(fitk,function(m) which.min(BIC(m))))
	print(k)
}

## OOS var
ll <- log(lambda)
cvm <- apply(D,2,mean)
cvs <- apply(D,2,sd)/sqrt(K-1)
cvlo <- cvm-cvs
cvhi <- cvm+cvs
seg.min <- which.min(cvm)
cv1se <- (cvm[seg.min]+cvs[seg.min])-cvm
seg.1se <- min((1:length(cvm))[cv1se>=0])


xdf <- rowMeans(sapply(fit,
	function(f) c(f$df,rep(tail(f$df,1),100-length(f$df)))))

library(glmnet)
fld <- c(1,rep.int(1:K,times=diff(chunks)))
fld[randi] <- fld
cvnet <- cv.glmnet(as.matrix(fgl[,1:9]), fgl$type, 
		family="multinomial",foldid=fld)

ylim <- c(1.8,3.1)#range(c(cvlo,cvhi,cvnet$cvm+c(-1,1)*cvnet$cvs),na.rm=TRUE)
## deviance 
pdf(file="fgl_cv.pdf",width=8,height=3.5)
par(mfrow=c(1,2),mai=c(1,.8,.8,.2),omi=c(0,.2,0,0))
plot(ll,cvm,type="n",xlab="",ylab="",
	ylim=ylim,bty="n")
mtext("dmr",side=3, font=2,line=2.5)
segments(x0=ll, y0=cvlo, y1=cvhi, col="grey70")
points(ll,cvm,pch=20,col=4)
abline(v=log(B@lambda[1]),col="darkorange",lty=2)
abline(v=log(Bbic@lambda[1]),col="green",lty=2)
abline(v=log(lambda[seg.min]), lty=2)
abline(v=log(lambda[seg.1se]), lty=2)
dfi <- unique(round(seq(1,100,length=ceiling(length(axTicks(1))))))
axis(3,at=log(lambda[dfi]), labels=round(xdf[dfi])-1,tick=FALSE)
par(xpd=FALSE)
lines(x=c(ll[100],ll[1]), y=rep(mean(icD[,'AIC']),2), lty=2, col="darkorange")
lines(x=c(ll[100],ll[1]), y=rep(mean(icD[,'BIC']),2), lty=2, col="green")
text(x=rep(ll[1],2),y=colMeans(icD),
	labels=c("aic","bic"),pos=4)

plot(cvnet,ylim=ylim,bty="n",xlab="", ylab="")
mtext("glmnet",side=3, font=2,line=2.5)
mtext("log lambda",side=1,font=3,outer=TRUE,line=-2)
mtext("Multinomial Deviance",side=2,font=3,outer=TRUE,line=-1/2)
dev.off()

R2 <- 1-D/D[,1]
print(colMeans(R2))
max(colMeans(R2))













