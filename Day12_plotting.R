library(testit)
library(klaR)

source('Day12_simulation.R')

PlotResults <- function(t,fitness_matrix)
{
# 	print("PlotResults: testing parameters")
# 	genotype_frequencies_adults <- CreateGenotypeFrequencies(0.1,0.2,0.7)
# 	fitness_matrix <- CreateFitnessMatrix1(1996)
# 	n_generations <- 5
#   t <- RunSimulation(genotype_frequencies_adults,fitness_matrix,n_generations)
# 	t
# 	print("~PlotResults: testing parameters")

	initial_A <- t[1,"AA"] + 0.5 * t[1,"Aa"]
	function_index <- 0
	year <- 0
	
	if (identical(fitness_matrix,CreateFitnessMatrix1(1996))) { year <- 1996; function_index <- 1 }
	if (identical(fitness_matrix,CreateFitnessMatrix2(1996))) { year <- 1996; function_index <- 2 }
	if (identical(fitness_matrix,CreateFitnessMatrix3(1996))) { year <- 1996; function_index <- 3 }
	if (identical(fitness_matrix,CreateFitnessMatrix1(2001))) { year <- 2001; function_index <- 1 }
	if (identical(fitness_matrix,CreateFitnessMatrix2(2001))) { year <- 2001; function_index <- 2 }
	if (identical(fitness_matrix,CreateFitnessMatrix3(2001))) { year <- 2001; function_index <- 3 }
	assert("",function_index!=0)
	assert("",year!=0)
	genotype_to_phenotype_function <- GetGenotypeToPhenotypeFunction(function_index)
	phenotype_to_genenotype_function <- GetPhenotypeToGenotypeFunction(function_index)
	genotype_to_phenotype_description <- GetGenotypeToPhenotypeFunctionDescription(function_index)
	
	
# 	png(filename=paste("Day12_",function_index,"_",initial_A * 100,"_",year,"_alleles_in_time.png",sep=""))
# 	plot(
# 		as.matrix(t$pA_gametes),ylim=c(0,1),
# 		type="l",col="red",
# 		main=paste("Allele frequencies in time for ",genotype_to_phenotype_description, " (",year,")"),
# 		xlab="Time (generations)",ylab="Allele frequency"
# 	)
# 	lines(as.matrix(t$pa_gametes),col="blue")
# 	legend_x <- 0.6 * n_generations
# 	legend_y <- 1.0
# 	legend(
# 		legend_x,
# 		legend_y,
# 		c("A","a"),col=c("red","blue"),pch=c(16,16)
# 	)
# 	dev.off()

	png(filename=paste("Day12_",function_index,"_",initial_A * 100,"_",year,"_genotypes_in_time.png",sep=""))
	plot(
		as.matrix(t$AA),ylim=c(0,1),
		type="l",col=PhenotypeToPlotColor(genotype_to_phenotype_function("AA")),
		main=paste("Genotype frequencies in time for ",genotype_to_phenotype_description, " (",year,")"),
		xlab="Time (generations)",ylab="Genotype frequency"
	)
	lines(as.matrix(t$Aa),col=PhenotypeToPlotColor(genotype_to_phenotype_function("Aa")))
	lines(as.matrix(t$aa),col=PhenotypeToPlotColor(genotype_to_phenotype_function("aa")))
	legend_x <- 0.6 * nrow(t)
	legend_y <- 1.0
	legend(
		legend_x,
		legend_y,
		c(
			paste("AA (",genotype_to_phenotype_function("AA"),")"),
			paste("Aa (",genotype_to_phenotype_function("Aa"),")"),
			paste("aa (",genotype_to_phenotype_function("aa"),")")
		),
		col=c(
			PhenotypeToPlotColor(genotype_to_phenotype_function("AA")),
			PhenotypeToPlotColor(genotype_to_phenotype_function("Aa")),
			PhenotypeToPlotColor(genotype_to_phenotype_function("aa"))
		),pch=c(16,16)
	)
  dev.off()

	png(filename=paste("Day12_",function_index,"_",initial_A * 100,"_",year,"_triplot.png",sep=""))
	triplot(as.matrix(t[,c("AA","Aa","aa")]),
		main=paste("Phenotypes in time for ",genotype_to_phenotype_description, " (",year,")"),
	  label=c(genotype_to_phenotype_function("AA"),genotype_to_phenotype_function("Aa"),genotype_to_phenotype_function("aa")),
	  grid = TRUE,
		type="l",
		col="red"
	)
	par(new=TRUE)
	triplot(as.matrix(CreatePredictedPhenotypes()),
		main=paste("Phenotypes in time for ",genotype_to_phenotype_description, " (",year,")"),
	  label=c(genotype_to_phenotype_function("AA"),genotype_to_phenotype_function("Aa"),genotype_to_phenotype_function("aa")),
	  grid = TRUE,
		type="l",
		col="black"
	)
  legend(
		0,
		0,
		c(
			"predicted",
			"observed"
		),
		col=c(
			"red",
			"black"
		),pch=c(16,16)
	)

  dev.off()

  return (t)
}

