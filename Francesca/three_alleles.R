w_sl<-matrix(data=c(1,2,0.3,1,1,4,2,0.1,1), nrow=3,ncol=3)  # 1996 payoff matrix
w_s<-matrix(data=c(1,0.7,2.2,2,1,0.9,1.2,1.8,1), nrow=3,ncol=3, byrow=TRUE) # 2001 payoff matrix

py<-0.2; pb<-0.3; po<-1-py-pb   # Probability value at time zero
p<-c(py,pb,po)
generations<-20
hob<-0.2
hby<-0.3
hyo<-0.1

sinervo_data<-matrix(data=c(0.326,0.568,0.106,
                            0.126,0.750,0.124,
                            0.361,0.328,  0.311,
                            0.420,0.321,  0.259,
                            0.485,0.389,  0.126,
                            0.350,0.551,  0.099,
                            0.426,0.452,	0.122,
                            0.603,0.150,	0.247,
                            0.524,0.374,	0.102,
                            0.452,0.460,	0.088),
                     ,nrow=10, ncol=3, byrow=TRUE)


# ITERATION EQUATION

three_alleles<-function(p,w_given,generations){
  p_alleles<-c(3*p[1]+p[2], 3*p[2]+p[3], 3*p[3]+p[1])/4
  w_phenotype<-w_given%*%p
  print(w_phenotype)
  w_alleles<-c((1-p_alleles[2])*w_phenotype[1]+p_alleles[2]*w_phenotype[2],
               (1-p_alleles[3])*w_phenotype[2]+p_alleles[3]*w_phenotype[3],
               (1-p_alleles[1])*w_phenotype[3]+p_alleles[1]*w_phenotype[1])

  w_bar<-sum(p_alleles*w_alleles)
  results<-matrix(0, nrow=generations, ncol=6)
  results[1,]<-c(p,p_alleles)
  
  for(i in 2:generations)
  {
    p_alleles<-(p_alleles*w_alleles)/w_bar
    w_alleles<-w_alleles*p_alleles/w_bar
    w_phenotype<-w_given%*%p
    w_bar<-sum(p_alleles*w_alleles)
    p<-c(p_alleles[1]^2+ 2*p_alleles[1]*p_alleles[2],
         p_alleles[2]^2+ 2*p_alleles[2]*p_alleles[3],
         p_alleles[3]^2+ 2*p_alleles[3]*p_alleles[1])
    results[i,]<-c(p,p_alleles)
  }
  return(results)
}

results_sl<-three_alleles(p,w_sl,generations)  # CONFRONT USING FIRST MATRIX 1996
results_s<-three_alleles(p,w_s,generations) # CONFRONT USING SECOND MATRIX 2001

library(klaR)  # Library to plot simplex
par(mfrow=c(1,2))
triplot(results_sl[,c(3,2,1)],label=c("O","B","Y"),pch=20,col=4, main="1996 fitness")
trilines(results_sl[,c(3,2,1)], col=4)
#tripoints(sinervo_data[,c(3,2,1)], pch=20,col=2)
#trilines(sinervo_data[,c(3,2,1)], col=2)

triplot(results_s[,c(3,2,1)],label=c("O","B","Y"),pch=20,col=4, main="2001 fitness")
trilines(results_s[,c(3,2,1)], col=4)
#tripoints(sinervo_data[,c(3,2,1)],pch=20,col=2)
#trilines(sinervo_data[,c(3,2,1)], col=2)

triplot(results_sl[,c(6,5,4)],label=c("O","B","Y"),pch=20,col=4, main="1996 fitness")
trilines(results_sl[,c(6,5,4)], col=4)
triplot(results_s[,c(6,5,4)],label=c("O","B","Y"),pch=20,col=4, main="2001 fitness")
trilines(results_s[,c(6,5,4)], col=4)