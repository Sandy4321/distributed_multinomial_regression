
library(distrom)
library(MASS)
library(glmnet)

data(fgl)
fit <- dmr(NULL, fgl[,1:9], fgl$type, 
	lambda.min.ratio=1e-4, gamma=0, tol=1e-8)

## plot the conditional MLE mu (which we've pegged at zero)
muhat <- log(1/rowSums(exp(predict(fit,fgl[,1:9]))))
pdf(file="fgl_mu.pdf", width=4,height=3)
par(mai=c(.9,.9,.2,.2))
hist(muhat,main="",xlab="Condtional MLE for mu",col="grey80")
dev.off()

## plot the individual Poisson model fit and selection
pdf(file="fgl_coef.pdf",width=8,height=4)
par(mfrow=c(2,3),mai=c(.3,.3,.6,.2),omi=c(.4,.4,0,0))
for(j in 1:6){
	plot(fit[[j]],xlab="",ylab="",col="grey70")
	mtext(names(fit)[j],font=2,line=2,cex=.9)
}
mtext("log lambda", outer=TRUE, font=3,side=1,line=1)
mtext("coefficient",outer=TRUE, font=3,side=2,line=1)
dev.off()

##### OOS experiment
K <- 20
n <- nrow(fgl)
chunks <- round(seq.int(1,n,length=K+1))
randi <- sample.int(n)
oosdev <- as.data.frame(matrix(0,ncol=5,nrow=K))
names(oosdev) <- c("net.CV1se","net.CVmin","dmr.CV1se","dmr.CVmin","dmr.AICc")

getdev <- function(f, lo, ...){ # multinomial deviance
			eta <- drop(predict(f,as.matrix(fgl[lo,1:9]), ...))
			ylo <- cbind(1:length(lo),fgl$type[lo])
			d <- eta[ylo] - log(rowSums(exp(eta))) 
			return(mean(-2*d)) }

for(k in 1:K){ # CV loop
	lo <- randi[chunks[k]:chunks[k+1]]
	dmrfit <- dmr(NULL, 
				fgl[-lo,1:9], fgl$type[-lo], lambda.min.ratio=1e-4)
	cvdmrfit <- dmr(NULL, cv=TRUE,
				fgl[-lo,1:9], fgl$type[-lo], lambda.min.ratio=1e-4)
	cvnet <- cv.glmnet(as.matrix(fgl[-lo,1:9]), fgl$type[-lo], family="multinomial")

	oosdev[k,"net.CV1se"] <- getdev(cvnet,lo,select="1se")
	oosdev[k,"net.CVmin"] <- getdev(cvnet,lo,select="min")
	oosdev[k,"dmr.AICc"] <- getdev(dmrfit,lo)
	oosdev[k,"dmr.CV1se"] <- getdev(cvdmrfit,lo,select="1se")
	oosdev[k,"dmr.CVmin"] <- getdev(cvdmrfit,lo,select="min")	

	print(k)
}

pdf(file="fgl_cv.pdf", width=7,height=3)
par(mai=c(.5,.9,.2,.2))
boxplot(oosdev, bty="n",col=c(rep("pink",2),rep("lightblue",3)),
	ylab="multinomial deviance",labels=c("net.CV1se"))
dev.off()







