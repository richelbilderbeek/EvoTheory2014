psi <- 0.5
omega <- 0.1
w_star <- function(z)
{
	p_star <- exp(-((z - psi)^2)/(2 * (omega ^ 2)))
} 

plot(w_star,xlim=c(0,2),main=paste("Optimal value: ",psi))
