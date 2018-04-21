mktable <-
function(cdata,ddata,rt="beta",varname,LG=1,Pv=0.00000005,Pc=0.979,Pd=0.979){

try(
 if(is.null(cdata))
  stop("No course data")
  )

try(
 if(is.null(ddata))
  stop("No disease data")
  )

cdata<-as.data.frame(cdata)
ddtat<-as.data.frame(ddata)

pvj<-cdata$pvj
pcj<-cdata$pcj

#print(length(posit))
newcdat<-subset(cdata,(pvj<=Pv&pcj>=Pc))
rsid<-newcdat$rsid
SNP.d<-ddata$SNP.d

idx<-is.element(SNP.d,rsid)
SNP<-SNP.d[idx]
ND<-ddata$N.case[idx]+ddata$N.ctr[idx]
pdj<-ND/max(ND)
newcdat<-as.data.frame(newcdat)

ddata<-cbind(ddata[idx,],pdj)
newcd<-fmerge(fl1=newcdat,fl2=ddata,ID1="rsid",ID2="SNP.d",A="",B="",method="no")
sbnewdat<-subset(newcd,pdj>=Pd)
sbnewdat<-as.data.frame(sbnewdat)

if(LG>2){

sbnewdat1<-sbnewdat[with(sbnewdat,order(chrn,posit)),]

chrn<-sbnewdat1$chrn
chrn.vect<-unique(chrn)
#ch<-chrn.vect[1]
posit<-sbnewdat1$posit

#distance<-posit[max(ix)]-posit[min(ix)]
distance<-rep(NA,length(chrn.vect))
posit.new<-rep(NA,length(posit))
for(ch in chrn.vect){
 #  print(ch)
    ix<-which(chrn==ch)
  #  print(ix)
    distance[ch]<-round((posit[max(ix)]-posit[min(ix)])/LG,0)
    if(is.na(distance[ch])==FALSE){
    if(distance[ch]<=1){
     posit.new[min(ix)]<-posit[min(ix)]
       }else{

      for(i in min(ix):(max(ix)-2)){
      if((posit[(i+1)]-posit[i]>LG)&(posit[(i+2)]-posit[(i+1)])>LG){
       posit.new[i]<-posit[i]
       posit.new[(i+1)]<-posit[(i+1)]
          }
        }
       }
    }
    }
# check if posit.new is in posit, T for yes and F for no
 indx<-is.element(posit,posit.new)
 sbnewdat2<-sbnewdat1
 sbnewdat3<-sbnewdat2[indx,]
 }else{
 sbnewdat3<-sbnewdat

 }

 dim(sbnewdat3)
sbnewdat3<-data.frame(sbnewdat3)
chr<-sbnewdat3$chrn
nSNP<-sbnewdat3$rsid
nposit<-sbnewdat3$posit

freqcase<-sbnewdat3$freq.case
 #print(is.numeric(freqcase))
 ln<-length(varname)
 rn<-nrow(sbnewdat3)
 betay<-alle<-freq<-sd<-matrix(NA,nrow=rn,ncol=ln)
 betay[,1]<-sbnewdat3$beta.d # disease beta value vector of SNPs
 alle[,1]<-sbnewdat3$a1.d
 freq[,1]<-sbnewdat3$freq.d
 sd[,1]<-sqrt(freqcase*(1-freqcase))
 for(i in 2:ln){
 k<-(2+2*(ln-1)+i)
 betay[,i]<-sbnewdat3[,k]
     k1<-(2+(ln-1)+i)
 freq[,i]<-sbnewdat3[,k1]
     k2<-2+i
 alle[,i]<-sbnewdat3[,k2]
 k3<-(2+3*(ln-1)+i)
 sd[,i]<-sbnewdat3[,k3]
 }
 #print(head(freq))
 rn<-nrow(betay)
 for(i in 2:ln){
 for(j in 1:rn){
 if(identical(alle[j,1],alle[j,i])==FALSE){
     betay[j,i]<--betay[j,i]
     }
   }
 }

colnames(betay)<-varname

betax<--matrix(NA,nrow=rn,ncol=ln)
idx<-seq_len(ln)
betax[,idx]<-betay[,idx]*sqrt((freq[,idx]*(1-freq[,idx]))/sd[,idx])

 betay<-cbind(nposit,betay)
 betay<-cbind(chr,nSNP,as.data.frame(betay))

 colnames(betay)[1:3]<-c("chr","SNP","posit")
 colnames(betax)<-varname
 betax<-cbind(chr,nSNP,nposit,betax)
 colnames(betax)[1:3]<-c("chr","SNP","posit")

rt<-tolower(rt)
if (rt=="beta"|rt=="b"){
 return(betay)
 }else if (rt=="path"|rt=="p"){
   return(betax)
 }

 }
