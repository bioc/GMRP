chrp <-
function(hg){
hg<-as.character(hg)	
HG<-unlist(strsplit(hg,split=":"))
n<-length(hg)
i<-seq(n)
k1<-i*2-1
k2<-i*2
chr<-HG[k1]
posit<-as.numeric(HG[k2])
CHR<-unlist(strsplit(chr,split="chr"))
chrn<-as.numeric(CHR[k2])
return(cbind(chrn,posit))
}
