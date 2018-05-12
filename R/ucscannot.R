ucscannot <-
function(UCSCannot,SNPn,A=3,B=1.9,C=1.3, method=1){
UCSCannot<-as.data.frame(UCSCannot)	
function_unit<-UCSCannot$function_unit
if(is.null(function_unit)){
stop("No function_unit found in the data")
}
genes<-UCSCannot$Symbol
if(is.null(genes)){
stop("No Symbol found in the data")
}
 N<-length(genes)
 gene.n<-length(unique(genes))

 intron<-subset(UCSCannot,function_unit=="intronic"|function_unit=="non-coding intronic")
 code<-subset(UCSCannot,function_unit=="coding")
 UTR3<-subset(UCSCannot,function_unit=="3utr")
 UTR5<-subset(UCSCannot,function_unit=="5utr")
 ncoding<-subset(UCSCannot,function_unit=="non-coding")
 upstream5<-subset(UCSCannot,function_unit=="5upstream")
 downstream3<-subset(UCSCannot,function_unit=="3downstream")

 intron<-length(intron[,1])/N
 coding<-length(code[,1])/N
 UTR3<-length(UTR3[,1])/N
 UTR5<-length(UTR5[,1])/N
 intergene<-length(ncoding[,1])/N
 upstream<-length(upstream5[,1])/N
 downstream<-length(downstream3[,1])/N

 res<-matrix(NA,1,8)
 res[1,]<-c(gene.n,coding,intron,UTR3,UTR5,intergene,upstream,downstream)
 colnames(res)<-c("genes","exons","introns","TUR3","UTR5","intergenes","upstream","downstream")
 res2<-c(coding,intron,intergene,UTR5,upstream,UTR3,downstream)
 oldcexmain<-par(cex.main=A)
 cols=c("gold2","red","magenta","blue","cornflowerblue","limegreen","mediumseagreen")
 main=paste("Distribution of",SNPn,"SNPs within",gene.n,"genes")
 
 if(method==1){
 annot<-c("exon","intron","intergene","5'UTR","5'upstream","3'UTR","3'downstream")
 lbls <- paste(annot, "\n", round(res2*100,1),"%", sep="")
# oldcexmain<-par(cex.main=A)
bisectors<-pie3D(res2,explode=0.1,main=main)
 pie3D.labels(radialpos=bisectors,radius=1,height=0.1,theta=pi/6,
       labels=lbls,labelcol=par("fg"),labelcex=B,labelrad=C,minsep=0.3)
  }else if(method==2){
        annot<-c("intron","exon","3'downstream","3'UTR","5'upstream","5'UTR","intergene")
          lbls <- paste(round(res2*100,1),"%", sep="")      
           # cols=c("gold2","red","magenta","blue","cornflowerblue","limegreen","mediumseagreen")
            bisectors<-pie3D(res2,explode=0.08,main=main)
            pie3D.labels(radialpos=bisectors,radius=1,height=0.1,theta=pi/6,
                  labels=lbls,labelcol=par("fg"),labelcex=(B+0.5),labelrad=C,minsep=0.3) 

if (B>=1.8){
 D<-1.8
 }else{
 D<-B}
 legend(0.5,1.05,annot,col=cols,pch=19,cex=D)
 }
 par(oldcexmain)
 return(res)
 }
