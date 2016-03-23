path <-
function(betav,model,outcome){
betav<-as.data.frame(betav)
rn<-nrow(betav) # number of beta values in each beta variable
cn<-ncol(betav)
#betav<-DataFrame(betav)
fit<-summary(lm(model,data=betav))

print(fit)
corr<-cor(betav)
cvar<-cov(betav)

beta.xx<-rep(0,(cn-1))

for(i in 2:cn){
beta.xx[i-1]<-fit$coef[i,1]
}

varx<-rep(0,(cn-1))
vary<-cvar[1,1]
corrx<-corr[-1,-1]
cvarx<-cvar[-1,-1]

for(i in 1:(cn-1)){
varx[i]<-cvarx[i,i]
}

py<-rep(0,(cn-1)) # direction path coefficients

for(i in 1:(cn-1)){
py[i]<-beta.xx[i]*(varx[i]/vary)^0.5
#f2[i]<-(varx[i]/vary)^0.5
}
print(py)
#print(vratio)
path<-matrix(0,(cn-1),(cn-1))
for(i in 1:(cn-1)){
	for(j in 1:(cn-1)){
path[i,j]<-corrx[i,j]*py[j]
 }
}
pe<-fit$sigma/sqrt(vary)

colnames(path)<-colnames(corrx)
rownames(path)<-rownames(corrx)

corr.outcome<-apply(path,1,sum)
path<-cbind(path,corr.outcome)
pn<-ncol(path)
beta<-fit$coef
se<-beta[,2]
#betacol<-matrix(NA,nrow=nrow(beta),ncol=1)
betar<-round(beta,6)
betar[,ncol(beta)]<-beta[,ncol(beta)] #p-values without round treatment
#print(betar)
corrname<-rownames(corr)
colnames(path)[pn]<-outcome

pathname<-colnames(path)

corr<-round(corr,6)
path<-round(path,6)
#R2<-rep(NA,cn)
#R2[1]<-fit$r.square
R_square<-rep(NA,ncol(path))
Pe<-rep(NA,ncol(path))
Pe[1]<-pe
#path<-rbind(path,R2)
R_square[1]<-round(fit$r.square,6)
Path_value<-rep(NA,length(py)+1)
Path_value[1:length(py)]<-round(py,6)
t_value<-p_value<-rep(NA,(1+length(py)))
for(i in 1:length(py)){
t_value[i]<-py[i]/se[i]
}
for(i in 1:length(py)){
p_value[i]<-2*pt(-abs(t_value[i]),df=rn-cn)
}

betax<-matrix(NA,nrow=nrow(betar),ncol=(cn-ncol(beta)))
betar<-cbind(betar,betax)

result<-rbind(betar,corrname,corr,pathname,path,Path_value,Pe,R_square)
return(noquote(result))

}
