pathdiagram <-
function(pathdata,disease,R2,range){

#library(diagram)

pathcad<-pathdata

#disease<-"CAD"
#range<-c(2:5)
pathname<-toupper(colnames(pathcad)[range])

rownames(pathcad)<-pathname

corr.xx<-pathcad[,range]
path_value<-pathcad$path
residual<-sqrt(1-R2)

par(mar = c(3, 1, 1, 1))
openplotmat()
shdv<-0.005
colorset<-c(colors()[645],colors()[642],colors()[503],colors()[373],colors()[92],colors()[134],colors()[139],colors()[39],colors()[20],colors()[44],colors()[71],colors()[123],colors()[128],colors()[132])
cl<-length(colorset)

mycolor<-rep(NA,cl)
mycolor<-colorset[1:cl]
mypath<-cbind(seq(1:length(path_value)),path_value)
mypath<-mypath[order(-path_value),]
indx<-mypath[,1]
mcolor<-rep(NA,(1+length(path_value)))

pl<-length(path_value)
for(i in 1:pl){
 if(mypath[i,2]<0){
 k=cl+1-i
 }else if(mypath[i,2]>0&mypath[i,2]<0.15){
 k=cl-5-i
 }else{
k=i
 }
mcolor[(indx[i]+1)]<-mycolor[k]
}
mcolor[1]<-colors()[82]

# Get plot coordinates
elpos <- coordinates(c(2, length(path_value)))
# adjust coordinates for Residual
elpos[2,1] <- abs((elpos[1,1]+elpos[1,2])/2)

# Specify Arrow positions
#1 Residual to Dependent
ft1 <- matrix(ncol = 2, byrow = TRUE, data = c(1, 2))
#2 Independent to dependent
ft2 <- matrix(ncol=2, byrow = FALSE,
data= c(seq((2+length(path_value)))[3:(length(path_value)+2)], rep(2, length(path_value))))

#3 For cor.x
fromto_CU <- t(combn(seq((2+length(path_value)))[3:(length(path_value)+2)],2))
#4 For path distances
fromto_ST <- rbind(ft1,ft2)

# Plot Path distance arrows
nr <- nrow(fromto_ST)

arrpos <- matrix(ncol = 2, nrow = nr)

for (i in 1:nr){
k=i+6
  arrpos[i, ] <-straightarrow (to = elpos[fromto_ST[i, 2], ],
                                from = elpos[fromto_ST[i, 1], ],
                     lcol=mcolor[i], lwd = 3, arr.pos = 0.7, arr.length = 0.5)
 }

#Label residual path distance arrow
text(arrpos[1, 1]-0.05, arrpos[1, 2] + 0.03,
     paste("Pe", " = ", round(residual, 3), sep=""), cex=1)

#Label path distance arrows
nr <- nrow(arrpos)
for(i in 2:nr){
  text(arrpos[i, 1]+(i*0.7-2)*0.14, arrpos[i, 2]-0.18,paste("P","(",pathname[i-1],")"," = ", round(path_value[i-1], 3), sep=""),cex=1,col=mcolor[i])

}

# Plot correlation arrows direction 1
nr <- nrow(fromto_CU)
arrpos <- matrix(ncol = 2, nrow = nr)
for (i in 1:nr){
#k=i+6
  arrpos[i, ] <- curvedarrow (to = elpos[fromto_CU[i, 2], ],
                                from = elpos[fromto_CU[i, 1], ],
                                lwd = 2, arr.pos = 0.8, arr.length = 0.5, curve = 0.3)
}

# Plot correlation arrows - direction 2
nr <- nrow(fromto_CU)
arrpos <- matrix(ncol = 2, nrow = nr)
for (i in 1:nr){
#k=i+6
  arrpos[i, ] <- curvedarrow (to = elpos[fromto_CU[i, 1], ],
                              from = elpos[fromto_CU[i, 2], ],
                              lwd = 2, arr.pos = 0.8, arr.length = 0.5, curve = -0.3)
}
# Create combinations of cor.x for labelling rxy in correlation arrows
rcomb <- as.data.frame(t(combn(seq(nrow(corr.xx)),2)))
rcomb <- paste(rcomb$V1,rcomb$V2, sep="")

# Label correlation arrows
nr <- nrow(fromto_CU)
arrpos <- matrix(ncol = 2, nrow = nr)
for (i in 1:nr)
#   k=i+6
  arrpos[i, ] <- curvedarrow (to = elpos[fromto_CU[i, 1], ],
                              from = elpos[fromto_CU[i, 2], ],
                              lwd = 2, arr.pos = 0.5, lcol = "transparent", arr.length = 0.5, curve = -0.3)

nr <- nrow(arrpos)
for(i in 1:nr){
  text(arrpos[i, 1], arrpos[i, 2] + 0.03,
       paste("r", rcomb[i]," = ", round(as.dist(corr.xx)[i], 3), sep=""), cex=1)
}

# Label Residual
#textrect (elpos[1,], 0.09, 0.03,lab = "Residual", box.col =colors()[(75+7)],
#          shadow.col = "grey", shadow.size = 0.01, cex = 1)

textround(elpos[1,], 0.02, 0.04,lab = "Residual", box.col =colors()[82],
          shadow.col = "grey", shadow.size = shdv, cex = 1)

# Label Dependent
#textrect (elpos[2,], 0.09, 0.03,lab = "AUDPC", box.col = "red",
#          shadow.col = "grey", shadow.size = 0.005, cex = 1)

textround (elpos[2,], 0.02, 0.04,lab = disease, box.col = "red",
          shadow.col = "grey", shadow.size = shdv, cex = 1)

# Label independents
nr <- nrow(elpos)
for (i in 3:nr){
	k<-i+5
 # textrect (elpos[i,], 0.09, 0.03,lab = colnames(x)[i-2], box.col =colors()[(75+k)],shadow.col = "grey", #shadow.size = 0.005, cex = 1)

 textround (elpos[i,], 0.02, 0.04,lab = pathname[i-2], box.col =mcolor[i-1],shadow.col = "grey", shadow.size = shdv, cex = 1)
}
}
