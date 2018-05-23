fmerge <-
function(fl1,fl2,ID1,ID2,A="",B="",method="No"){
fl1<-DataFrame(fl1)
fl2<-DataFrame(fl2)
try (if(is.null(dim(fl1))|is.null(dim(fl2))) 
   stop("No data to be merged")
   )
fl1<-cbind(fl1[,ID1],as.data.frame(fl1))
colnames(fl1)[1]<-"SNPID"
fl2<-cbind(fl2[,ID2],as.data.frame(fl2))
colnames(fl2)[1]<-"SNPID"
method<-tolower(method)
if(method=="no"|method=="n"){
merge.dat<-merge(fl1,fl2,by="SNPID",suffixes=c(A,B),sort=TRUE)
#merge.dat<-subset(merge.dat,SNPID!=".")
}
else if(method=="file1"){
merge.dat<-merge(fl1,fl2,by="SNPID",suffixes=c(A,B),all.x=TRUE,sort=TRUE)
}
else if(method=="file2"){
merge.dat<-merge(fl1,fl2,by="SNPID",suffixes=c(A,B),all.y=TRUE,sort=TRUE)
}
else if(method=="all"|method=="a"){
colnames(fl1)[1]<-"SNPID1"
colnames(fl2)[2]<-"SNPID2"
rownames(fl1)<-fl1$SNPID1
rownames(fl2)<-fl2$SNPID2
pos.all <- unique(c(fl1$SNPID1,fl2$SNPID2))
ind.fl1 <- match(pos.all,fl1$SNPID1)
ind.fl2 <- match(pos.all,fl2$SNPID2)
fl1.new <- fl1[ind.fl1,]
fl2.new<-fl2[ind.fl1,]
colnames(fl1.new) <- paste(A,colnames(fl1),sep='.')
colnames(fl2.new) <- paste(B,colnames(fl2),sep='.')
merge.dat<-cbind(as.data.frame(fl1.new),as.data.frame(fl2.new))
mergdata<-merge.dat[,-"SNPID2"]
}
mergdata<-merge.dat[,-1]
mergdata<-DataFrame(mergdata)
return(mergdata)
}
