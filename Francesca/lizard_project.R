# INITALIZE STARTING VALUES

w_sl<-matrix(data=c(1,1,2,2,1,0.1,0.3,4,1), nrow=3,ncol=3, byrow=TRUE)  # 1996 payoff matrix
w_s<-matrix(data=c(1,0.7,2.2,2,1,0.9,1.2,1.8,1), nrow=3,ncol=3, byrow=TRUE) # 2001 payoff matrix

py<-0.6; pb<-0.2; po<-1-py-pb   # Probability value at time zero
p<-c(py,pb,po)
generations<-100

# REPLICATOR EQUATION FUNCTION

replicator<-function(p,w,w_sl,generations){
  w<-w_sl%*%p
  w_bar<-sum(p*w)
  results<-matrix(0, nrow=generations, ncol=7)
  
  for(i in 1:generations)
  {
      p<-(p*w)/w_bar
      w<-w_sl%*%p
      w_bar<-sum(p*w)
      results[i,]<-c(p,w,w_bar)
  }
return(results)
}

sinervo_data<-matrix(data=c(0.326,0.568,0.106,
                            0.126,0.750,0.124,
                            0.361,0.328,  0.311,
                            0.420,0.321,  0.259,
                            0.485,0.389,	0.126,
                            0.350,0.551,	0.099,
                            0.426,0.452,	0.122,
                            0.603,0.150,	0.247,
                            0.524,0.374,	0.102,
                            0.452,0.460,	0.088),
                     ,nrow=10, ncol=3, byrow=TRUE)


results_sl<-replicator(p,w,w_sl,generations)  # CONFRONT USING FIRST MATRIX 1996
results_s<-replicator(p,w,w_s,generations) # CONFRONT USING SECOND MATRIX 2001

par(mfrow=c(1,2))

# FITNESS PLOT IN TIME
plot(results_sl[,1], type="l", col=6, ylim=c(0,1), xlab="generations", ylab="frequencies")
lines(results_sl[,2], col=4)
lines(results_sl[,3], col=2)
plot(results_s[,1], type="l", col=6, ylim=c(0,1), xlab="generations", ylab="frequencies")
lines(results_s[,2], col=4)
lines(results_s[,3], col=2)

# FITNESS SIMPLEX PLOT

library(klaR)  # Library to plot simplex
triplot(results_sl[,c(3,2,1)],label=c("O","B","Y"),pch=20,col=4, main="1996 fitness")
trilines(results_sl[,c(3,2,1)], col=4)
tripoints(sinervo_data[,c(3,2,1)], pch=20,col=2)
trilines(sinervo_data[,c(3,2,1)], col=2)


triplot(results_s[,c(3,2,1)],label=c("O","B","Y"),pch=20,col=4, main="2001 fitness")
trilines(results_s[,c(3,2,1)], col=4)
tripoints(sinervo_data[,c(3,2,1)],pch=20,col=2)
trilines(sinervo_data[,c(3,2,1)], col=2)
