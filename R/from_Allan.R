
BIN.RAST<-matrix(nrow=22,ncol=17)

for(i in 1:17)
{
  for(j in 1:22)
  {
    TIB.IJ<-TIB.BeechJulyJoined %>% filter(PSIclass==i,Quantile==j)
    BIN.RAST[23-j,i]<-nrow(TIB.IJ)
  }
  print(i)
}

BIN.RAST.APP<-apply(BIN.RAST,2,FUN=function(x){x/sum(x)})

RAST<-rast(BIN.RAST.APP[,2:16])

ext(RAST)<-c(-16,-1,1,22)

COL<-colorRampPalette(c("grey","yellow","red"))

plot(RAST,col=COL(100))

##and here the paper

###https://onlinelibrary.wiley.com/doi/10.1111/gcb.15057