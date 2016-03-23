snpposit <-
function(SNPdata,SNP_hg19="chr",X="no",LG= 10,main="A",maxd=2000){
SNPdata<-as.data.frame(SNPdata)	
if(SNP_hg19!="chr"){
hg19<-as.character(SNPdata$SNP_hg19)
#print(hg19[1:10])

chromp<-chrp(hg19)
#print(chromp[1:6,])
}else{
chr<-SNPdata$chr
if(length(chr)==0){
stop("no data")
}

SNPdata<-SNPdata[with(SNPdata,order(chr,posit)),]
chrn<-as.numeric(SNPdata$chr)
posit<-as.numeric(SNPdata$posit)
chromp<-cbind(chrn,posit)

#print(CHR)
}
chr<-chromp[,1]
CHR<-unique(chr)
nch<-length(CHR) # number of chromosomes
meanl<-m<-md<-rep(NA,nch)

#calculate SNP number on each chromosome
nsnpchr<-function(x,LG){
y<-as.numeric(x[,2])
#print(y)
nr<-length(y)
y<-sort(y)
m<-0
ix<-seq(nr-1)	
iy<-seq(from=2,to=nr)
m<-length(which(y[iy]-y[ix]>=LG))
return(m)
}

#averaged distance between SNPs on each chromosome
mdist<-function(x,LG){
x<-as.data.frame(x)
y<-x[,2]
#print(y)
nr<-length(y)
 y<-sort(y)
 mdist<-(y[nr]-y[1])/((nr-1)*LG)
return(mdist)

}

#calculate average distance between SNPs and SNP number on each chromosome
m<-md<-snpn<-rep(0,nch)
for(i in 1:nch){
chrmp<-subset(chromp,chrn==i)
snpn[i]<-nrow(chrmp)
if(snpn[i]!=0){
m[i]<-nsnpchr(subset(chromp,chrn==i),LG)
md[i]<-mdist(subset(chromp,chrn==i),LG)
  }
}

md[which(md>=maxd)]<-maxd
md[which(is.na(md))]<-0
#print(length(md))
#print(length(CHR))
colors<-c("red","lightsalmon","lightseagreen","mediumspringgreen","mediumorchid4","mediumpurple4","lightskyblue4","mediumseagreen","deeppink",
"mediumvioletred","firebrick","magenta","navyblue","lightsteelblue2","navajowhite","blue","maroon4","orange2","purple","darkorchid1",
"black","cyan")

chrcol<-rep(NA,nch)
for(i in 1:nch){
 chrcol[i]<-colors[i]
}
if(X=="yes"){
CHR[nch]<-"X"	
}
#print(CHR)
par(oma=c(2,1,2,1))
bp<-barplot(md,width=2,ylab="mean length (kbp) of interval between SNPs", names.arg =paste("chr",CHR,sep=""),col=chrcol,ylim=c(0,(maxd+250)),
main=main,cex.lab=1.5,cex.names=1,cex.main=2.5,cex.axis=1.5)

# space between lab and column
YY<-rep(NA,nch)
for(i in 1:nch){
	if(!is.na(md[i])){
if(md[i]>=1000){
YY[i]<-md[i]*1.05
   }else if(md[i]<1000&md[i]>=500){
   YY[i]<-md[i]*1.08
   }else if(md[i]<500&md[i]>=200){
   YY[i]<-md[i]*1.2
   }else if(md[i]<200&md[i]>=100){
    YY[i]<-(md[i]*1.8)
    }else if(md[i]<100&md[i]>=50){
    YY[i]<-(md[i]*2.5)
    }else if(md[i]<50&md[i]>0){
    YY[i]<-(md[i]*5)
    }else if(md[i]==0){
    YY[i]<-45
  }
  }
}
#text(x=bp,y=YY,labels=n,col="black")
text(x=bp,y=YY,labels=snpn,col="black",cex=1.3)
return(m)
}
