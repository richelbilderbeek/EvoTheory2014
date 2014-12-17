rm(list = ls())
options(error = browser)
options(echo = FALSE)
library(testit)

# source('SinervoData.R')
source('Day12_test.R')
source('Day12_plotting.R')
source('Day12_simulation.R')

TestDay12()

run_all_sims <- FALSE
if (run_all_sims == TRUE)
{
	for (year in c(1996,2001))
	{
		for (function_index in seq(1:3))
		{
  		for (initial_phenotype_densities in c(CreatePhenotypeFrequencies(0.1,0.6),CreatePredictedPhenotypesRaw()[1]))
	  	{
			  n_generations <- 50
	      if (year == 2001) { n_generations <- 100 }
		    t <- RunSimulation(initial_phenotype_densities,function_index,n_generations,year)
	      PlotResults(t,initial_phenotype_densities,function_index,n_generations,year)
		  }
		}
	}
}

TestSimulation()
# t <- RunSimulation(CreatePhenotypeFrequencies(0.1,0.6),1,50,1996)
# t
# PlotResults(t,0.05,1,50,1996)
