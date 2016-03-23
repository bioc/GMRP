pathdiagram2 <-
function(pathD,pathO,rangeD,rangeO,disease,R2D,R2O){
#require(agricolae)
#require(Hmisc)
#library(diagram)

pathcad<-pathD
pathtc<-pathO
#rangeD<-c(2:5)
#rangeO<-c(2:4)
#R2D<-0.241792
#R2O<-0.965284
cadpathname<-toupper(colnames(pathcad)[rangeD])

rownames(pathcad)<-cadpathname
rownametc<-toupper(colnames(pathtc)[rangeO])
for(i in 1:length(rownametc))
if(identical(cadpathname[i],rownametc[i])==FALSE){
	stop("no match between child and parent causal variables.")
}
if(R2D==""|R2O==""){
stop("R2D and R2O are required")
}
pathtc<-data.frame(pathtc)
pathcad<-data.frame(pathD)
#rownames(pathtc)<-rownametc
corr.tc<-pathtc[,rangeO]
colnames(corr.tc)<-rownametc
rownames(corr.tc)<-rownametc
path.tc<-pathtc$path
path.cad<-pathcad$path
residual.tc<-sqrt(1-R2O)
residual.cad<-sqrt(1-R2D)

par(mar = c(1, 1, 1, 1))
openplotmat()
# Get plot coordinates
elpos <- coordinates(c(2, length(path.tc)))*0.7

# adjust coordinates for Residual
elpos[2,1] <- abs((elpos[1,1]+elpos[1,2])/2)

# Specify Arrow positions
#1 Residual to Dependent
ft1 <- matrix(ncol = 2, byrow = TRUE, data = c(1, 2))
#2 Independent to dependent
ft2 <- matrix(ncol=2, byrow = FALSE,
data= c(seq((2+length(path.tc)))[3:(length(path.tc)+2)], rep(2, length(path.tc))))
#3 For cor.x
fromto_CU <- t(combn(seq((2+length(path.tc)))[3:(length(path.tc)+2)],2))
#4 For path distances
fromto_ST <- rbind(ft1,ft2)

# Plot Path distance arrows
nr <- nrow(fromto_ST)

arrpos <- matrix(ncol = 2, nrow = nr)

for (i in 1:nr){
k=i+6
  arrpos[i, ] <-straightarrow (to = elpos[fromto_ST[i, 2], ],
                                from = elpos[fromto_ST[i, 1], ],
                     lcol=colors()[(75+k)], lwd = 2, arr.pos = 0.6, arr.length = 0.5)
 }

#Label residual path distance arrow
text(arrpos[1, 1]-0.05, arrpos[1, 2] + 0.02,
     paste("Pe", " = ", round(residual.tc, 3), sep=""), cex=0.8)

#Label path distance arrows
nr <- nrow(arrpos)
for(i in 2:nr){
 k=i+6
text(arrpos[i, 1]+(i-2)*0.06, arrpos[i, 2]-0.05,paste("P","(",tolower(rownametc[i-1]),")"," = ", round(path.tc[i-1], 3), sep=""),cex=0.8,col=colors()[(75+k)])

}

# Plot correlation arrows direction 1
nr <- nrow(fromto_CU)
arrpos <- matrix(ncol = 2, nrow = nr)
for (i in 1:nr){
#k=i+6
  arrpos[i, ] <- curvedarrow (to = elpos[fromto_CU[i, 2], ],
                                from = elpos[fromto_CU[i, 1], ],
                                lwd = 2, arr.pos = 0.8, arr.length = 0.5, curve = 0.35)
}

# Plot correlation arrows - direction 2
for (i in 1:nr){
#k=i+6
  arrpos[i, ] <- curvedarrow (to = elpos[fromto_CU[i, 1], ],
                                from = elpos[fromto_CU[i, 2], ],
                                lwd = 2, arr.pos = 0.8, arr.length = 0.5, curve = -0.35)
}

# Create combinations of cor.x for labelling rxy in correlation arrows
rcomb <- as.data.frame(t(combn(seq(nrow(corr.tc)),2)))
rcomb <- paste(rcomb$V1,rcomb$V2, sep="")

# Label correlation arrows
nr <- nrow(fromto_CU)
arrpos <- matrix(ncol = 2, nrow = nr)
for (i in 1:nr)
#   k=i+6
  arrpos[i, ] <- curvedarrow (to = elpos[fromto_CU[i, 1], ],
                              from = elpos[fromto_CU[i, 2], ],
                              lwd = 2, arr.pos = 0.5, lcol = "transparent", arr.length = 0.5, curve = -0.35)

nr <- nrow(arrpos)
for(i in 1:nr){
  text(arrpos[i, 1], arrpos[i, 2] + 0.03,
       paste("r", rcomb[i]," = ", round(as.dist(corr.tc)[i], 2), sep=""), cex=1)
}

# Label Residual
#textrect (elpos[1,], 0.09, 0.03,lab = "Residual", box.col =,
#          shadow.col = "grey", shadow.size = 0.005, cex = 1)

textround(c(elpos[1,1]-0.08,elpos[1,2]), 0.02, 0.04,lab = "Residual", box.col =colors()[(75+7)],
          shadow.col = "grey", shadow.size = 0.005, cex = 1)

# Label Dependent
#textrect (elpos[2,], 0.09, 0.03,lab = "AUDPC", box.col = "red",
#          shadow.col = "grey", shadow.size = 0.005, cex = 1)

# the second outcome variable boxround
textround (elpos[2,], 0.02, 0.04,lab = cadpathname[length(cadpathname)], box.col = colors()[122],
          shadow.col = "grey", shadow.size = 0.005, cex = 1.3)

# Label independents
nr <- nrow(elpos)
for (i in 3:nr){
	k<-i+5
 # textrect (elpos[i,], 0.09, 0.03,lab = colnames(x)[i-2], box.col =colors()[(75+k)],shadow.col = "grey", #shadow.size = 0.005, cex = 1)

 textround (elpos[i,], 0.02, 0.04,lab = colnames(corr.tc)[i-2], box.col =colors()[(75+k)],shadow.col = "grey", shadow.size = 0.005, cex = 1)
}

elpos1 <- coordinates(c(2, length(path.cad)))*1.15

# adjust coordinates for Residual
elpos1[2,1] <- abs((elpos1[1,1]+elpos1[1,2])*0.75)

# Specify Arrow positions
#1 Residual to Dependent
fit1 <- matrix(ncol = 2, byrow = TRUE, data = c(1, 2))
#2 Independent to dependent
fit2 <- matrix(ncol=2, byrow = FALSE,
data= c(seq((2+length(path.cad)))[3:(length(path.cad)+2)], rep(2, length(path.cad))))
#3 For cor.x
fromto_CU1 <- t(combn(seq((2+length(path.cad)))[3:(length(path.cad)+2)],2))
#4 For path distances
fromto_ST1 <- rbind(fit1,fit2)

# Plot Path distance arrows
nr1 <- nrow(fromto_ST1)-1

arrpos1 <- matrix(ncol = 2, nrow = nr1)

arrpos1[1, ] <-straightarrow (to = c(elpos1[fromto_ST1[1, 2], 1]-0.5,elpos1[fromto_ST1[1, 2], 2]),
                                from = c(elpos1[fromto_ST1[1, 1], 1]-0.1,elpos1[fromto_ST1[1, 1], 2]),
                     lcol=colors()[(75+7)], lwd = 2, arr.pos = 0.5, arr.length = 0.5)

 # path line of the second outcome (TC) to the fist outcome variable (disease)
arrpos1[2, ] <-straightarrow (to = c(elpos1[fromto_ST1[1, 2], 1]-0.5,elpos1[fromto_ST1[1, 2], 2]),
                                from = c(elpos[fromto_ST[1, 2], 1],elpos[fromto_ST[1, 2], 2]+0.03),
                     lcol=colors()[122], lwd = 2, arr.pos = 0.6, arr.length = 0.7)

for (i in 2:nr1){
if(i==3){
k=i+6
  arrpos1[i, ] <-curvedarrow(to = c(elpos1[fromto_ST1[i, 2], 1]-0.5,elpos1[fromto_ST1[i, 2], 2]),
                                from = elpos[fromto_ST[i, 1], ],
                     lcol=colors()[(75+k)], lwd = 2, arr.pos = 0.7, arr.length = 0.5, curve = 0.05)
}else{
k=i+6
  arrpos1[i, ] <-straightarrow (to = c(elpos1[fromto_ST1[i, 2], 1]-0.5,elpos1[fromto_ST1[i, 2], 2]),
                                from = elpos[fromto_ST[i, 1], ],
                     lcol=colors()[(75+k)], lwd = 2, arr.pos = 0.7, arr.length = 0.5)
 }
}
#Label residual path distance arrow

text(arrpos1[1, 1]-0.04, arrpos1[1, 2] + 0.03,
     paste("Pe", " = ", round(residual.cad, 3), sep=""), cex=0.8)

# residual coefficient to the firest outcome(disease)
textround(c(elpos1[1,1]-0.18,elpos1[1,2]), 0.02, 0.04,lab = "Residual", box.col =colors()[(75+7)],
          shadow.col = "grey", shadow.size = 0.005, cex = 1)

textround (c(elpos1[2,1]-0.5,elpos1[2,2]), 0.02, 0.04,lab = disease, box.col = "red",
          shadow.col = "grey", shadow.size = 0.005, cex = 1)

 # path coefficient of second outcome variable (such as TC)to first outcome variable(disease)
 text(arrpos1[2, 1]+0.1, arrpos1[2, 2]+0.03,
     paste("P","(", cadpathname[length(cadpathname)],")", " = ", round(path.cad[length(path.cad)], 3), sep=""), cex=0.8,col=colors()[122])

nr1 <- nrow(arrpos1)
for(i in 2:nr1){
 k=i+6
text(arrpos1[i, 1]+(i*0.8-2)*0.09, arrpos1[i, 2]-(i-1)*0.05,paste("P","(", cadpathname[i-1],")"," = ", round(path.cad[i-1], 3), sep=""),cex=1,col=colors()[(75+k)])

}

# recover the color and labs
nr <- nrow(elpos)
for (i in 3:nr){
	k<-i+5
 # textrect (elpos[i,], 0.09, 0.03,lab = colnames(x)[i-2], box.col =colors()[(75+k)],shadow.col = "grey", #shadow.size = 0.005, cex = 1)

 textround (elpos[i,], 0.02, 0.04,lab = colnames(corr.tc)[i-2], box.col =colors()[(75+k)],shadow.col = "grey", shadow.size = 0.005, cex = 1.3)
}

}
