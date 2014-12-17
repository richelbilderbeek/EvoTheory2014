library(testit)
library(klaR)

PlotResults <- function(t,initial_A,function_index,n_generations,year)
{
	genotype_to_phenotype_function <- GetGenotypeToPhenotypeFunction(function_index)
	phenotype_to_genenotype_function <- GetPhenotypeToGenotypeFunction(function_index)
	genotype_to_phenotype_description <- GetGenotypeToPhenotypeFunctionDescription(function_index)

	png(filename=paste("Day12_",function_index,"_",initial_A * 100,"_",year,"_alleles_in_time.png",sep=""))
	plot(
		as.matrix(t$pA_gametes),ylim=c(0,1),
		type="l",col="red",
		main=paste("Allele frequencies in time for ",genotype_to_phenotype_description, " (",year,")"),
		xlab="Time (generations)",ylab="Allele frequency"
	)
	lines(as.matrix(t$pa_gametes),col="blue")
	legend_x <- 0.6 * n_generations
	legend_y <- 1.0
	legend(
		legend_x,
		legend_y,
		c("A","a"),col=c("red","blue"),pch=c(16,16)
	)
	dev.off()

	png(filename=paste("Day12_",function_index,"_",initial_A * 100,"_",year,"_genotypes_in_time.png",sep=""))
	plot(
		as.matrix(t$pAA_juveniles),ylim=c(0,1),
		type="l",col=PhenotypeToPlotColor(genotype_to_phenotype_function("AA")),
		main=paste("Genotype frequencies in time for ",genotype_to_phenotype_description, " (",year,")"),
		xlab="Time (generations)",ylab="Genotype frequency"
	)
	lines(as.matrix(t$pAa_juveniles),col=PhenotypeToPlotColor(genotype_to_phenotype_function("Aa")))
	lines(as.matrix(t$paa_juveniles),col=PhenotypeToPlotColor(genotype_to_phenotype_function("aa")))
	legend_x <- 0.6 * n_generations
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
	triplot(as.matrix(t[,c("pY","pB","pO")]),
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