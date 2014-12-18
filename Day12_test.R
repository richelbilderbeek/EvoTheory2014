library(testit)

source('Day12_simulation.R')
source('Day12_plotting.R')

TestDay12 <- function()
{
	# Basic tests
#   {
#   	a <- matrix(data = c(-1,-2,-3,-4),2,2)
#   	b <- matrix(data = c(1,2,3,4),2,2)
#   	c <- matrix(data = c(11,22,33,44),2,2)
#   	t <- c(a,b,c)
#   	assert("I would not expect this to fail",t[1] == a)
#   }
  {
  	m <- data.frame(data = matrix(0,3,3))
  	names(m) = c("aa","Aa","AA")
  	rownames(m) = c("Aa","aa","AA")
    m <- SortFitnessMatrix(m)
  	assert("",names(m) == c("AA","Aa","aa"))
  	assert("",rownames(m) == c("AA","Aa","aa"))
  }
	{
	  t <- CreateFitnessMatrix1(1996)
		# t[rare,common]
	  assert("",t["AA","Aa"]==0.3)
	  assert("",t["AA","aa"]==4.0)
	  assert("",t["aa","AA"]==0.1)
	  assert("",t["Aa","AA"]==2.0)
	}

  {
	  t <- CreateFitnessMatrix2(1996)
		# t[rare,common]
	  assert("",t["Aa","aa"]==0.3)
	  assert("",t["Aa","AA"]==4.0)
	  assert("",t["AA","Aa"]==0.1)
	  assert("",t["aa","Aa"]==2.0)
	}

  {
	  t <- CreateFitnessMatrix3(1996)
		# t[rare,common]
	  assert("",t["aa","AA"]==0.3)
	  assert("",t["aa","Aa"]==4.0)
	  assert("",t["Aa","aa"]==0.1)
	  assert("",t["AA","aa"]==2.0)
	}
	#	CalculateNextGenotypeDensities
  {
  	genotype_frequencies_parents <- CreateGenotypeFrequencies(0.1,0.2,0.7)
  	genotype_frequencies_parents
  	fitness_matrix <- CreateFitnessMatrix1(1996)
	  genotype_frequencies_offspring <- CalculateNextGenotypeDensities(genotype_frequencies_parents,fitness_matrix)
  	genotype_frequencies_offspring
  	assert("",genotype_frequencies_offspring != genotype_frequencies_parents)
  }
  #	RunSimulation
  {
  	genotype_frequencies_adults <- CreateGenotypeFrequencies(0.1,0.2,0.7)
  	fitness_matrix <- CreateFitnessMatrix1(1996)
  	n_generations <- 5
	  t <- RunSimulation(genotype_frequencies_adults,fitness_matrix,n_generations)
  	t
  }
	# PlotResults
  {
  	genotype_frequencies_adults <- CreateGenotypeFrequencies(0.1,0.2,0.7)
  	fitness_matrix <- CreateFitnessMatrix1(1996)
  	n_generations <- 5
	  t <- RunSimulation(genotype_frequencies_adults,fitness_matrix,n_generations)
		PlotResults(t,fitness_matrix)
  }
}

TestDay12()