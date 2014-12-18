rm(list = ls())
options(error = browser)
options(echo = FALSE)
library(testit)

source('Day12_test.R')
source('Day12_plotting.R')
source('Day12_simulation.R')

TestDay12()

DoIt <- function()
{
	for (year in c(1996,2001))
	{
		for (function_index in c(1,2,3))
		{
			function_index <- 1
			year <- 1996
		  fitness_matrix <- 0
			if (function_index == 1) fitness_matrix <- CreateFitnessMatrix1(year)
			if (function_index == 2) fitness_matrix <- CreateFitnessMatrix2(year)
			if (function_index == 3) fitness_matrix <- CreateFitnessMatrix3(year)	
  	  assert("",ncol(fitness_matrix) == nrow(fitness_matrix))

			for (initial_genotype_frequencies_adults_index in c(1,2))
			{
				initial_genotype_frequencies_adults <- 0
				if (initial_genotype_frequencies_adults_index == 1) initial_genotype_frequencies_adults <- CreateGenotypeFrequencies(0.1,0.6)
				if (initial_genotype_frequencies_adults_index == 2) initial_genotype_frequencies_adults <- CreatePredictedPhenotypesRaw()[1,]
				
				assert("",ncol(initial_genotype_frequencies_adults) == nrow(fitness_matrix))
	
				n_generations <- 50
			  t <- RunSimulation(initial_genotype_frequencies_adults,fitness_matrix,n_generations)
				PlotResults(t,fitness_matrix)
			}
		}
	}
}

DoIt()
