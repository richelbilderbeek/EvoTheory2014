rm(list = ls())
options(error = browser)
options(echo = FALSE)

library(testit)
library(tableplot)
library(rgl)
library(colorRamps)
library(lattice)

fitnesses <- data.frame(0.0,3,3)
fitnesses[1,1] <- 0.808
fitnesses[1,2] <- 1.000
fitnesses[1,3] <- 0.842
fitnesses[2,1] <- 0.974	
fitnesses[2,2] <- 0.997	
fitnesses[2,3] <- 0.636
fitnesses[3,1] <- 0.916	
fitnesses[3,2] <- 0.682	
fitnesses[3,3] <- 0.393
rownames(fitnesses) <- c("BB","Bb","bb")
colnames(fitnesses) <- c("AA","Aa","aa")
fitnesses
inbreeding_coefficient <- 0.0
linkage_disequilibrium <- 0.0

# image.plot(fitnesses)

assert("inbreeding_coefficient must be e [0.0,1.0]",inbreeding_coefficient >= 0.0 && inbreeding_coefficient < 1.0)
assert("linkage_disequilibrium must be e [-0.25,0.25]",linkage_disequilibrium >= -0.25 && linkage_disequilibrium < 0.25)
# Show the fitness matrix in orange, because orange is a color with high fitness
# persp3d(data.matrix(fitnesses),col = "orange")

n_fitnesses <- 100
fitness_landscape <- data.frame(
	matrix(0.0,n_fitnesses,n_fitnesses),
	row.names=seq(0,1,1/(n_fitnesses-1))
)
names(fitness_landscape) = seq(0,1,1/(n_fitnesses-1))
for (x in seq(1,n_fitnesses))
{	
	pA <- (x-1) / (n_fitnesses - 1)
	assert("pA must be e [0.0,1.0]",pA >= 0.0 && pA <= 1.0)
	pa <- 1.0 - pA

	# Calculate the genotype frequencies for the A/a locus
	genotype_frequencies_a <- data.frame(rep(0,3),row.names=c("AA","Aa","aa"))
	names(genotype_frequencies_a) = c("p")
  genotype_frequencies_a[1,1] <- ((1.0-F) * (pA * pA)) + (F * pA)
  genotype_frequencies_a[2,1] <- ((1.0-F) * (2.0 * pA * pa))
  genotype_frequencies_a[3,1] <- ((1.0-F) * (pa * pa)) + (F * pa)
	assert("genotype_frequencies_a must add up to one",abs(sum(genotype_frequencies_a) - 1.0) < 0.001)
  for (y in seq(1,n_fitnesses))
	{
		pB <- (y-1) / (n_fitnesses - 1)
		assert("pB must be e [0.0,1.0]",pB >= 0.0 && pB <= 1.0)
		pb <- 1.0 - pB
    # Calculate the genotype frequencies for the B/b locus
  	genotype_frequencies_b <- data.frame(rep(0,3),row.names=c("BB","Bb","bb"))
	  names(genotype_frequencies_b) = c("p")
	  genotype_frequencies_b[1,1] <- ((1.0-F) * (pB * pB)) + (F * pB)
	  genotype_frequencies_b[2,1] <- ((1.0-F) * (2.0 * pB * pb))
	  genotype_frequencies_b[3,1] <- ((1.0-F) * (pb * pb)) + (F * pb)
		assert("genotype_frequencies_b must add up to one",abs(sum(genotype_frequencies_b) - 1.0) < 0.001)
		
		#Calculate the fitness values
		w <-
			( genotype_frequencies_a[1,1]  
		  * (
		  	    genotype_frequencies_b[1,1] * fitnesses[1,1]
		  	  + genotype_frequencies_b[2,1] * fitnesses[2,1]
		  	  + genotype_frequencies_b[3,1] * fitnesses[3,1]
		    )
			)
		  +	
			( genotype_frequencies_a[2,1]  
		  * (
		  	    genotype_frequencies_b[1,1] * fitnesses[1,2]
		  	  + genotype_frequencies_b[2,1] * fitnesses[2,2]
		  	  + genotype_frequencies_b[3,1] * fitnesses[3,2]
		    )
			)
		  +	
			( genotype_frequencies_a[3,1]  
		  * (
		  	    genotype_frequencies_b[1,1] * fitnesses[1,3]
		  	  + genotype_frequencies_b[2,1] * fitnesses[2,3]
		  	  + genotype_frequencies_b[3,1] * fitnesses[3,3]
		    )
			)
    fitness_landscape[y,x] <- w
	}
}

fitness_table <- data.frame(
	matrix(0.0,nrow(fitness_landscape)*ncol(fitness_landscape),3)
)
names(fitness_table) <- c("x","y","w")

for (y in c(1:nrow(fitness_landscape)))
{
  for (x in c(1:ncol(fitness_landscape)))
  {
    i <- 1 + ( (y-1) * ncol(fitness_landscape) ) + (x-1)
  	assert("",i >= 1 && i <= nrow(fitness_landscape)*ncol(fitness_landscape))
  	fitness_table[i,1] <- names(fitness_landscape)[x]
  	fitness_table[i,2] <- row.names(fitness_landscape)[y] 
  	fitness_table[i,3] <- fitness_landscape[y,x]
  }
}
	

#persp3d(data.matrix(fitness_landscape),col = "orange")
#splot(data.matrix(fitness_landscape),col = "orange")
#c = cut(seq(0,1,0.1), breaks=10)
# cols = blue2red(100)[as.numeric(c)]
plot3d(data.matrix(fitness_table),xlab="f(A)",ylab="f(B)",zlab="Fitness")
levelplot(fitness_table[,3]~fitness_table[,1]*fitness_table[,2] ,
	data=fitness_table,contour=TRUE,
	xlab="f(A)",ylab="f(B)",zlab="Fitness"
)
#splot(data.matrix(fitness_landscape),col = "orange")

combinations<-cbind(combinations,meanfitness)
levelplot(meanfitness~combinations[,1]*combinations[,2],data=combinations,contour=TRUE,col.regions=blue2green2red(5000),xlab="Frequency of A",ylab="Frequency of C")



print("DONE")




stop()
stop()
stop()
stop()
stop()

# CalculateAlleleFitnesses <- function(allele_frequencies,fitnesses,inbreeding_coefficient = 0.0)
# {
# 	assert("Gene frequencies must add up to one", abs(sum(allele_frequencies) - 1.0) < 0.0001)
# 	assert("Inbreeding_coefficient must be e [0,1]",inbreeding_coefficient >= 0.0);
# 	assert("Inbreeding_coefficient must be e [0,1]",inbreeding_coefficient <= 1.0);
#   assert("Might improve this later", length(allele_frequencies) == 2)
#   assert("Might improve this later", length(fitnesses) == 3)
#   p <- allele_frequencies[1]
#   q <- allele_frequencies[2]
# 	ic <- inbreeding_coefficient
#   Ws <- fitnesses
# 	allele_fitnesses <- c(
# 		(ic * Ws[1]) + ((1-ic) * p * Ws[1]) + ((1-ic) * q * Ws[2]),
# 		(ic * Ws[3]) + ((1-ic) * q * Ws[3]) + ((1-ic) * p * Ws[2])
#   )
# 	return (allele_fitnesses)
# }
# 
# CalculateMeanFitnesses <- function(allele_frequencies,allele_fitnesses)
# {
#   assert("There must be as much allele frequencies as allele fitnesses",
#   	length(allele_frequencies) == length(allele_fitnesses))
#   mean_fitness = 	(allele_frequencies %*% allele_fitnesses)[1][1]
# 	# assert("",length(mean_fitness) == 1)
# 	return (mean_fitness)
# }
# 
# 
# CreateGeneFrequencies <- function(allele_frequencies,inbreeding_coefficient = 0.0)
# {
# 	assert("CreateGeneFrequencies: allele_frequencies must be a vector",is.vector(allele_frequencies))
# 	if (abs(sum(allele_frequencies) - 1.0) > 0.0001) { browser() }
# 	assert("CreateGeneFrequencies: Allele frequencies must add up to one", abs(sum(allele_frequencies) - 1.0) < 0.0001)
# 	assert("CreateGeneFrequencies: Inbreeding_coefficient must be e [0,1]",inbreeding_coefficient >= 0.0);
# 	assert("CreateGeneFrequencies: Inbreeding_coefficient must be e [0,1]",inbreeding_coefficient <= 1.0);
#   assert("CreateGeneFrequencies: Might improve this later", length(allele_frequencies) == 2)
# 	p <- allele_frequencies[1]
# 	q <- allele_frequencies[2]
# 	ic <- inbreeding_coefficient
# 	gene_frequencies <- c(
# 	  ((1.0 - ic) * p * p      ) + (ic * p),
# 		((1.0 - ic) * p * q * 2.0)          ,
# 		((1.0 - ic) * q * q      ) + (ic * q)
#   );
# 	return (gene_frequencies)
# }
# 
# # Gene frequencies
# assert("Two allele frequencies must yield three gene frequences", length(CreateGeneFrequencies(c(0,1))) == 3);
# assert("When only one allele is present, all gene frequencies must be the homozygote", CreateGeneFrequencies(c(1,0)) == c(1,0,0));
# assert("When only one allele is present, all gene frequencies must be the homozygote", CreateGeneFrequencies(c(0,1)) == c(0,0,1));
# assert("p = q = 0.5 should be p*p, 2*p*q, q*q",CreateGeneFrequencies(c(0.5,0.5)) == c(0.25,0.5,0.25));
# assert("p = 0.1 = 0.9 should be p*p, 2*p*q, q*q",CreateGeneFrequencies(c(0.1,0.9)) - c(0.01,0.18,0.81) < rep(0.0000001,3));
# assert("Gene frequencies must sum up to one",sum(CreateGeneFrequencies(c(0.0,1.0))) == 1);
# assert("Gene frequencies must sum up to one", sum(CreateGeneFrequencies(c(1.0,0.0))) == 1);
# assert("Gene frequencies must sum up to one", sum(CreateGeneFrequencies(c(0.1,0.9))) == 1);
# assert("At full inbreeding, there should be no heterozygotes",CreateGeneFrequencies(c(0.5,0.5),1.0) == c(0.5,0.0,0.5));
# assert("At full inbreeding, there should be no heterozygotes",CreateGeneFrequencies(c(0.1,0.9),1.0) == c(0.1,0.0,0.9));
# assert("At no inbreeding, Wp = p*Wpp + q*Wpq",CalculateAlleleFitnesses(c(0.3,0.7),c(0.7,0,0.3),0.0) == c(0.21,0.21))
# assert("At full inbreeding, Wp = Wpp and Wq = Wqq",CalculateAlleleFitnesses(c(0.5,0.5),c(2,0,3),1.0) == c(2.0,3.0))
# 
# CalculateNextAlleleFrequencies <- function(allele_frequencies,fitnesses,inbreeding_coefficient)
# {
# 	gene_frequencies <- CreateGeneFrequencies(allele_frequencies)
# 	allele_fitnesses <- CalculateAlleleFitnesses(allele_frequencies,fitnesses,inbreeding_coefficient)
# 	mean_fitness <- CalculateMeanFitnesses(allele_frequencies,allele_fitnesses)
# 	next_allele_frequencies <- allele_frequencies * allele_fitnesses / mean_fitness	
# 	return (next_allele_frequencies)
# }
# 
# 
# Exercise6 <-function()
# {
# 	fitnesses_males   <- c(2.0,0.5,1.0)
# 	fitnesses_females <- c(1.0,0.5,2.0)
# 	inbreeding_coefficient <- 0.0
# 	initial_allele_frequencies_males <- c(0.9,0.1)
# 	initial_allele_frequencies_females <- c(0.9,0.1)
# 	
# 	allele_frequencies <- matrix(0,0,length(initial_allele_frequencies_males) + length(initial_allele_frequencies_females))
# 	t <- c(initial_allele_frequencies_males,initial_allele_frequencies_female)
# 	allele_frequencies <- rbind(allele_frequencies,initial_allele_frequencies)
# 	for (i in c(1:10))
# 	{
# 	  allele_frequencies <- rbind(
# 	  	allele_frequencies,
# 	  	CalculateNextAlleleFrequencies(allele_frequencies[i,],fitnesses,inbreeding_coefficient)
# 	  )
# 	}
#   plot(allele_frequencies[,1],ylim=c(0,1),type="l",col="red",ylab="Allele frequency",xlab="Time",main="Exercise 5B")
# 	lines(allele_frequencies[,2],col="blue")
# }
# 
# 
# Exercise6()
# 
# browser()
# 
# 
# t <- c(1,2)
# u <- c(3,4)
# w <- c(t,u)
# w

