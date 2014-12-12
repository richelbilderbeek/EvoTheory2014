rm(list = ls())
options(error = browser)
library(testit)

CalculateAlleleFitnesses <- function(allele_frequencies,fitnesses,inbreeding_coefficient = 0.0)
{
	assert("Gene frequencies must add up to one", abs(sum(allele_frequencies) - 1.0) < 0.0001)
	assert("Inbreeding_coefficient must be e [0,1]",inbreeding_coefficient >= 0.0);
	assert("Inbreeding_coefficient must be e [0,1]",inbreeding_coefficient <= 1.0);
  assert("Might improve this later", length(allele_frequencies) == 2)
  assert("Might improve this later", length(fitnesses) == 3)
  p <- allele_frequencies[1]
  q <- allele_frequencies[2]
	ic <- inbreeding_coefficient
  Ws <- fitnesses
	allele_fitnesses <- c(
		(ic * Ws[1]) + ((1-ic) * p * Ws[1]) + ((1-ic) * q * Ws[2]),
		(ic * Ws[3]) + ((1-ic) * q * Ws[3]) + ((1-ic) * p * Ws[2])
  )
	return (allele_fitnesses)
}

CalculateMeanFitnesses <- function(allele_frequencies,allele_fitnesses)
{
  assert("There must be as much allele frequencies as allele fitnesses",
  	length(allele_frequencies) == length(allele_fitnesses))
  mean_fitness = 	(allele_frequencies %*% allele_fitnesses)[1][1]
	# assert("",length(mean_fitness) == 1)
	return (mean_fitness)
}


CreateGeneFrequencies <- function(allele_frequencies,inbreeding_coefficient = 0.0)
{
	assert("CreateGeneFrequencies: allele_frequencies must be a vector",is.vector(allele_frequencies))
	if (abs(sum(allele_frequencies) - 1.0) > 0.0001) { browser() }
	assert("CreateGeneFrequencies: Allele frequencies must add up to one", abs(sum(allele_frequencies) - 1.0) < 0.0001)
	assert("CreateGeneFrequencies: Inbreeding_coefficient must be e [0,1]",inbreeding_coefficient >= 0.0);
	assert("CreateGeneFrequencies: Inbreeding_coefficient must be e [0,1]",inbreeding_coefficient <= 1.0);
  assert("CreateGeneFrequencies: Might improve this later", length(allele_frequencies) == 2)
	p <- allele_frequencies[1]
	q <- allele_frequencies[2]
	ic <- inbreeding_coefficient
	gene_frequencies <- c(
	  ((1.0 - ic) * p * p      ) + (ic * p),
		((1.0 - ic) * p * q * 2.0)          ,
		((1.0 - ic) * q * q      ) + (ic * q)
  );
	return (gene_frequencies)
}

# Gene frequencies
assert("Two allele frequencies must yield three gene frequences", length(CreateGeneFrequencies(c(0,1))) == 3);
assert("When only one allele is present, all gene frequencies must be the homozygote", CreateGeneFrequencies(c(1,0)) == c(1,0,0));
assert("When only one allele is present, all gene frequencies must be the homozygote", CreateGeneFrequencies(c(0,1)) == c(0,0,1));
assert("p = q = 0.5 should be p*p, 2*p*q, q*q",CreateGeneFrequencies(c(0.5,0.5)) == c(0.25,0.5,0.25));
assert("p = 0.1 = 0.9 should be p*p, 2*p*q, q*q",CreateGeneFrequencies(c(0.1,0.9)) - c(0.01,0.18,0.81) < rep(0.0000001,3));
assert("Gene frequencies must sum up to one",sum(CreateGeneFrequencies(c(0.0,1.0))) == 1);
assert("Gene frequencies must sum up to one", sum(CreateGeneFrequencies(c(1.0,0.0))) == 1);
assert("Gene frequencies must sum up to one", sum(CreateGeneFrequencies(c(0.1,0.9))) == 1);
assert("At full inbreeding, there should be no heterozygotes",CreateGeneFrequencies(c(0.5,0.5),1.0) == c(0.5,0.0,0.5));
assert("At full inbreeding, there should be no heterozygotes",CreateGeneFrequencies(c(0.1,0.9),1.0) == c(0.1,0.0,0.9));
assert("At no inbreeding, Wp = p*Wpp + q*Wpq",CalculateAlleleFitnesses(c(0.3,0.7),c(0.7,0,0.3),0.0) == c(0.21,0.21))
assert("At full inbreeding, Wp = Wpp and Wq = Wqq",CalculateAlleleFitnesses(c(0.5,0.5),c(2,0,3),1.0) == c(2.0,3.0))

CalculateNextAlleleFrequencies <- function(allele_frequencies,fitnesses,inbreeding_coefficient)
{
	gene_frequencies <- CreateGeneFrequencies(allele_frequencies)
	allele_fitnesses <- CalculateAlleleFitnesses(allele_frequencies,fitnesses,inbreeding_coefficient)
	mean_fitness <- CalculateMeanFitnesses(allele_frequencies,allele_fitnesses)
	next_allele_frequencies <- allele_frequencies * allele_fitnesses / mean_fitness	
	return (next_allele_frequencies)
}


Exercise6 <-function()
{
	fitnesses_males   <- c(2.0,0.5,1.0)
	fitnesses_females <- c(1.0,0.5,2.0)
	inbreeding_coefficient <- 0.0
	initial_allele_frequencies_males <- c(0.9,0.1)
	initial_allele_frequencies_females <- c(0.9,0.1)
	
	allele_frequencies <- matrix(0,0,length(initial_allele_frequencies_males) + length(initial_allele_frequencies_females))
	t <- c(initial_allele_frequencies_males,initial_allele_frequencies_female)
	allele_frequencies <- rbind(allele_frequencies,initial_allele_frequencies)
	for (i in c(1:10))
	{
	  allele_frequencies <- rbind(
	  	allele_frequencies,
	  	CalculateNextAlleleFrequencies(allele_frequencies[i,],fitnesses,inbreeding_coefficient)
	  )
	}
  plot(allele_frequencies[,1],ylim=c(0,1),type="l",col="red",ylab="Allele frequency",xlab="Time",main="Exercise 5B")
	lines(allele_frequencies[,2],col="blue")
}


Exercise6()

browser()


t <- c(1,2)
u <- c(3,4)
w <- c(t,u)
w
