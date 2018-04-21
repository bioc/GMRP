snpPositAnnot<-function(SNPdata,SNP_hg19="chr",main){
#This function provide mean length of interval between SNP positions on each chromosome
#and number of SNPs on each chromosome and ouput numbers of SNPs with interval length>=LG.
#The data contain hg19 or chr and position. If SNP_hg19 == chr, then data must contain chr column and posit column
#if SNP_hg19 != chr, then SNP_hg19 must be hg19 or hg18 that contain a column of "chrx:xxxxxxx". 

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
SNPdata<-SNPdata[order(chr),]
#print(SNPdata[1:3,])
chrn<-as.numeric(SNPdata$chr)
posit<-as.numeric(SNPdata$posit)
chromp<-cbind(chrn,posit)

#print(CHR)
}
chr<-chromp[,1]
CHR<-unique(chr)
nch<-length(CHR) # number of chromosomes

#print(nch)

meanl<-m<-md<-rep(NA,nch)

mdist<-function(x){

y<-as.numeric(x[,2])
#print(y)
nr<-length(y)
y<-sort(y)
if(nr==0){
lgdist<-0	
}else if(nr==1){
lgdist<-log(y[1])		
	}else{
dist<-rep(NA,nr-1)
for(i in 1:nr-1){	
dist[i]<-(y[(i+1)]-y[i])
}
lgdist<-log(min(dist))
#print(dist)
}
return(lgdist)		
	
}

m<-md<-snpn<-rep(0,nch)
for(i in 1:nch){
chrmp<-subset(chromp,chrn==CHR[i])
snpn[i]<-nrow(chrmp)
if(snpn[i]!=0){
#	m[i]<-0
#	md[i]<-0}else{

#m[i]<-nsnpchr(subset(chromp,chrn==i))

md[i]<-mdist(subset(chromp,chrn==CHR[i]))

}

}


#mlength<-meanl[1:20]
colors<-c("red","lightsalmon","lightseagreen","mediumspringgreen","mediumorchid4","mediumpurple4","lightskyblue4","mediumseagreen","deeppink",
"mediumvioletred","firebrick","magenta","navyblue","lightsteelblue2","navajowhite","blue","maroon4","orange2","purple","darkorchid1",
"black","cyan")

chrcol<-rep(NA,nch)
for(i in 1:nch){
 chrcol[i]<-colors[i]	
}

par(oma=c(2,1,2,1))
bp<-barplot(md,width=2,ylab="log least length (kbp) of intervals between SNPs", names.arg =paste("chr",CHR,sep=""),col=chrcol,ylim=c(0,30),
main=main,cex.lab=2.0,cex.names=2,cex.main=2.5,cex.axis=1.5)


# space between lab and column
YY<-rep(NA,nch)
for(i in 1:nch){
	if(!is.na(md[i])){
if(md[i]>=10){	
  YY[i]<-md[i]+1.5
   }else if(md[i]<10&md[i]>0){
    YY[i]<-md[i]+1.5
    }else if(md[i]==0){
   YY[i]<-2	
   }
  }
}
#text(x=bp,y=YY,labels=n,col="black")
# for data SNP247
text(x=bp,y=YY,labels=snpn,col="black",cex=2.5)
# for data SNP368
#return(m)
}
