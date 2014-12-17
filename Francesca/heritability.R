results<-matrix(0, nrow=99, ncol=2)

# CYCLE CHANGING ALLELE FREQUENCY

for(i in 1:99){
  p1<-i*0.01  # FREQUENCY OF ALLELE 1
  p2<-1-p1
  p_alleles<-c(p1,p2)
  
  # MATRIX CONTAINING PROBABILITY OF THE 3X3 GRID
  
  mat<-matrix(data=c(0,p1*p2^2,p2^3,
                     p2*p1^2,p1*p2,p1^2*p2,
                     p1^3,p2*p1^2,0), nrow=3, ncol=3, byrow=TRUE)

  x<-c(-4,0,1)  # GENOTYPE COORDINATES. ARBITRARY!!
  y<-mat%*%x    # Y IS OFFSPRING COORDINATES
  h<-coef(lm(y ~ x))[2]   #FIT LINE
  results[i,]<-c(p1,h)
}
plot(results[,1],-2*results[,2], xlab="p_a", ylab="h^2", type="l", col=4)
