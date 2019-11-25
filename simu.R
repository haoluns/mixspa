load("first3_sparse.RData")


library(fda);library(dplyr)
library(Rmpfr)


library(dplyr)
library(reshape2)
library(data.table)
library(lsei)
source("funs_mixspa.R")


library(fda)
timepts=1:10;
norder=4 ## cubic B-spline
nbasis=norder+length(timepts)-2;
spline_basis=create.bspline.basis(rangeval=c(1,10),nbasis=6,norder=4)

library(fda)
library(fdapace)

soapfpca<-function(observedx,timepointsx){
res_pace<- FPCA(observedx, timepointsx,list(dataType='Sparse',error=TRUE, kernel='epan', 
	verbose=FALSE,methodBwCov="GCV",methodBwMu="GCV"))
coef_mat0 = coef(Data2fd(argvals = res_pace$workGrid,y=res_pace$phi,spline_basis))

pc1s = first_FPC(coef_mat0[,1],observed=observedx, timepoints=timepointsx,minit=6,gamma=1e1,threshold=1e-5)

previous_beta = list()
previous_beta[[1]] = pc1s$beta
pc2s = third_FPC_conditional2(rnorm(spline_basis$nbasis), observed=pc1s$residual, timepoints=timepointsx, gamma=3e4,betalist =previous_beta,threshold=1e-4,range_score=max(abs(pc1s$sfit)))

previous_beta[[2]] = pc2s$beta

pc3s = third_FPC_conditional2(rnorm(spline_basis$nbasis), observed=pc2s$residual, timepoints=timepointsx, gamma=1e2,betalist =previous_beta,threshold=1e-3,range_score=max(abs(pc2s$sfit)))

return(list(pc1=pc1s,pc2=pc2s,pc3=pc3s))
}



mixFPCA<-function(observed,timepoints){

	sigma = 10

	fmeanfit = findmean(observed, timepoints, minit=6,gamma=0,threshold=5e-3)

	observedmu1 = lapply(1:length(observed),function(i){
		 (eval.fd(timepoints[[i]], fmeanfit$pc_fit1))[,1]
	})
	observedmu2 = lapply(1:length(observed),function(i){
		 (eval.fd(timepoints[[i]], fmeanfit$pc_fit2))[,1]
	})

	
     

	if (mean(sapply(observedmu1,mean))<mean(sapply(observedmu2,mean))){
	tempp  = observedmu1
	observedmu1 = observedmu2
	observedmu2 = tempp
	
	mu1 = fmeanfit$pc_fit2
	mu2 = fmeanfit$pc_fit1
	fmeanfit_sigmak = rev(fmeanfit$sigmak)
	reverseflag = 1
	}else{
		mu1 = fmeanfit$pc_fit1
	mu2 = fmeanfit$pc_fit2
	fmeanfit_sigmak=fmeanfit$sigmak
	reverseflag = 0
	}
	
	
	observedcentered1 = lapply(1:length(observed),function(i){
		 observed[[i]]-(eval.fd(timepoints[[i]], mu1))[,1]
	})
	observedcentered2 = lapply(1:length(observed),function(i){
		 observed[[i]]-(eval.fd(timepoints[[i]], mu2))[,1]
	})


	owik=meanfit$wikmat
	if(reverseflag == 0){
	classest = sapply(owik[,2],function(x)ifelse(x>0.5,1,2))
	}else{
	classest = sapply(owik[,2],function(x)ifelse(x>0.5,2,1))
	}
	
	observed1 = observedcentered1[which(classest==1)]
	timepoints1 = timepoints[which(classest==1)]
	observed2 = observedcentered2[which(classest==2)]
	timepoints2 = timepoints[which(classest==2)]

	
	res_pace1<- FPCA(observed1, timepoints1,list(dataType='Sparse',error=TRUE, kernel='epan', verbose=TRUE,methodBwCov="GCV",methodBwMu="GCV",methodSelectK = 3))

	res_pace2<- FPCA(observed2, timepoints2,list(dataType='Sparse',error=TRUE, kernel='epan', verbose=TRUE,methodBwCov="GCV",methodBwMu="GCV",methodSelectK = 3))
	

	pacefpc11 = Data2fd(argvals = res_pace1$workGrid,y=res_pace1$phi[,1],spline_basis)
	pacefpc12 = Data2fd(argvals = res_pace1$workGrid,y=res_pace1$phi[,2],spline_basis)
	pacefpc13 = Data2fd(argvals = res_pace1$workGrid,y=res_pace1$phi[,3],spline_basis)
	
	pacefpc21 = Data2fd(argvals = res_pace2$workGrid,y=res_pace2$phi[,1],spline_basis)
	pacefpc22 = Data2fd(argvals = res_pace2$workGrid,y=res_pace2$phi[,2],spline_basis)
	pacefpc23 = Data2fd(argvals = res_pace2$workGrid,y=res_pace2$phi[,3],spline_basis)


	
   fpc1x = soapfpca(observed1,timepoints1)
   fpc2x = soapfpca(observed2,timepoints2)



	return(list(
	mu1 = mu1,
	mu2 = mu2,
	
	pc11 = fpc1x$pc1$pc_fit,
	pc12 = fpc2x$pc1$pc_fit,

	pc21 = fpc1x$pc2$pc_fit,
	pc22 = fpc2x$pc2$pc_fit,
	
	pc31 = fpc1x$pc3$pc_fit,
	pc32 = fpc2x$pc3$pc_fit,
	
	meanfit = fmeanfit,
	
	pacepc11 = pacefpc11,
	pacepc12 = pacefpc12,
	pacepc13 = pacefpc13,
	pacepc21 = pacefpc21,
	pacepc22 = pacefpc22,
	pacepc23 = pacefpc23	
	
	))

}


ssize = 150

timegrid = seq(1,10,by=0.05)
pp=c(0.3,0.4,0.5)
sparse=TRUE

for(jj in 1:3){
for(it in 1:100){


clusters = lapply(1:ssize,function(i){
x= rbinom(1,1,pp[jj])
ifelse(x==1,1,2)
})



timepoints = lapply(1:ssize,function(i){
curtimegrid = sample(timegrid,20)
sort(curtimegrid)
})


observed = lapply(1:ssize,function(i){
	if(sparse==TRUE){
	tm=timepoints[[i]]
	}else{
	tm = timegrid
	}
	mu1 = eval.fd(tm, meanfit$pc_fit1) 
	mu2 = eval.fd(tm, meanfit$pc_fit2)
    
	
	pc1add1 = eval.fd(tm, fpc1$pc1$pc_fit)*rnorm(1,0,50)
	pc1add2 = eval.fd(tm, fpc2$pc1$pc_fit)*rnorm(1,0,50)

	pc2add1 = eval.fd(tm, fpc1$pc2$pc_fit)*rnorm(1,0,25)
	pc2add2 = eval.fd(tm, fpc2$pc2$pc_fit)*rnorm(1,0,25)
	
	pc3add1 = eval.fd(tm, fpc1$pc3$pc_fit)*rnorm(1,0,10)
	pc3add2 = eval.fd(tm, fpc2$pc3$pc_fit)*rnorm(1,0,10)


	err1=rnorm(1,mean=0,sd=0.1)
	err2=rnorm(1,mean=0,sd=0.1)

	if(clusters[i]==1){
	return(as.numeric(mu1+pc1add1+pc2add1+pc3add1+err1))
	}else{
	return(as.numeric(mu2+pc1add2+pc2add2+pc3add1+err2))
	}

})

if(sparse==FALSE){
timepoints = lapply(1:ssize,function(i){
timegrid
})
}

#inprod(pacemu-meanfit$pc_fit,pacemu-meanfit$pc_fit)

res_cfpca = mixFPCA(observed,timepoints)
errorcfpca_mu1 = inprod(res_cfpca$mu1-meanfit$pc_fit1,res_cfpca$mu1-meanfit$pc_fit1)
errorcfpca_mu2 = inprod(res_cfpca$mu2-meanfit$pc_fit2,res_cfpca$mu2-meanfit$pc_fit2)

errorcfpca_pc11= min(inprod(res_cfpca$pc11-fpc1$pc1$pc_fit,res_cfpca$pc11-fpc1$pc1$pc_fit),inprod(res_cfpca$pc11+fpc1$pc1$pc_fit,res_cfpca$pc11+fpc1$pc1$pc_fit))
errorcfpca_pc12= min(inprod(res_cfpca$pc12-fpc2$pc1$pc_fit,res_cfpca$pc12-fpc2$pc1$pc_fit),inprod(res_cfpca$pc12+fpc2$pc1$pc_fit,res_cfpca$pc12+fpc2$pc1$pc_fit))

errorcfpca_pc21 = min(inprod(res_cfpca$pc21-fpc1$pc2$pc_fit,res_cfpca$pc21-fpc1$pc2$pc_fit),inprod(res_cfpca$pc21+fpc1$pc2$pc_fit,res_cfpca$pc21+fpc1$pc2$pc_fit))
errorcfpca_pc22 = min(inprod(res_cfpca$pc22-fpc2$pc2$pc_fit,res_cfpca$pc22-fpc2$pc2$pc_fit),inprod(res_cfpca$pc22+fpc2$pc2$pc_fit,res_cfpca$pc22+fpc2$pc2$pc_fit))


errorcfpca_pc31 = min(inprod(res_cfpca$pc31-fpc1$pc3$pc_fit,res_cfpca$pc31-fpc1$pc3$pc_fit),inprod(res_cfpca$pc31+fpc1$pc3$pc_fit,res_cfpca$pc31+fpc1$pc3$pc_fit))
errorcfpca_pc32 = min(inprod(res_cfpca$pc32-fpc2$pc3$pc_fit,res_cfpca$pc32-fpc2$pc3$pc_fit),inprod(res_cfpca$pc32+fpc2$pc3$pc_fit,res_cfpca$pc32+fpc2$pc3$pc_fit))




errorpace_pc11= min(inprod(res_cfpca$pacepc11-fpc1$pc1$pc_fit,res_cfpca$pacepc11-fpc1$pc1$pc_fit),inprod(res_cfpca$pacepc11+fpc1$pc1$pc_fit,res_cfpca$pacepc11+fpc1$pc1$pc_fit))
errorpace_pc12= min(inprod(res_cfpca$pacepc21-fpc2$pc1$pc_fit,res_cfpca$pacepc21-fpc2$pc1$pc_fit),inprod(res_cfpca$pacepc21+fpc2$pc1$pc_fit,res_cfpca$pacepc21+fpc2$pc1$pc_fit))

errorpace_pc21 = min(inprod(res_cfpca$pacepc12-fpc1$pc2$pc_fit,res_cfpca$pacepc12-fpc1$pc2$pc_fit),inprod(res_cfpca$pacepc12+fpc1$pc2$pc_fit,res_cfpca$pacepc12+fpc1$pc2$pc_fit))
errorpace_pc22 = min(inprod(res_cfpca$pacepc22-fpc2$pc2$pc_fit,res_cfpca$pacepc22-fpc2$pc2$pc_fit),inprod(res_cfpca$pacepc22+fpc2$pc2$pc_fit,res_cfpca$pacepc22+fpc2$pc2$pc_fit))


errorpace_pc31 = min(inprod(res_cfpca$pacepc13-fpc1$pc3$pc_fit,res_cfpca$pacepc13-fpc1$pc3$pc_fit),inprod(res_cfpca$pacepc13+fpc1$pc3$pc_fit,res_cfpca$pacepc13+fpc1$pc3$pc_fit))
errorpace_pc32 = min(inprod(res_cfpca$pacepc23-fpc2$pc3$pc_fit,res_cfpca$pacepc23-fpc2$pc3$pc_fit),inprod(res_cfpca$pacepc23+fpc2$pc3$pc_fit,res_cfpca$pacepc23+fpc2$pc3$pc_fit))




o_sigma = res_cfpca$meanfit$sigma
o_pik = res_cfpca$meanfit$pik
owik = res_cfpca$meanfit$wikmat

classest = ifelse(owik[,2]>0.5,1,2)

classtrue = clusters%>%do.call(c,.)

correctrate = max(length(which(classest==classtrue))/ssize,1-length(which(classest==classtrue))/ssize)


output = c(	
errorcfpca_mu1,errorcfpca_mu2,
errorcfpca_pc11,errorcfpca_pc12,
errorcfpca_pc21,errorcfpca_pc22,
errorcfpca_pc31,errorcfpca_pc32,
errorpace_pc11,errorpace_pc12,
errorpace_pc21,errorpace_pc22,
errorpace_pc31,errorpace_pc32,
o_sigma,o_pik,correctrate
)



if (it==1){
outputmat = output
}else{
outputmat = rbind(outputmat,output)
}

}


}

