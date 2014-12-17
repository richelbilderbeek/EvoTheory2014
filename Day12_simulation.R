library(testit)

source('SinervoData.R')
source('../R/Genetics.R')

CreatePhenotypeFrequencies <- function(pY, pB, pO = 1.0 - pY - pB)
{
	n_rows <- 1
	n_cols <- 3
	fs <- data.frame(matrix(0,n_rows,n_cols))
	fs[1,1] <- pY
	fs[1,2] <- pB
	fs[1,3] <- pO
	rownames(fs) <- c("f")
	colnames(fs) <- c("Y","B","O")
	assert("Phenotype Frequencies must sum up to one or be empty",
		abs(sum(fs) - 1.0) < 0.0001
	)
	return (fs)
}

GenotypeToPhenotype1 <- function(genotype)
{
	if (genotype == "AA") return ("O")
	if (genotype == "Aa") return ("Y")
	if (genotype == "aa") return ("B")
	assert("Unknown genotype, should be 'AA', 'Aa' or 'aa'",1==2)
}

GenotypeToPhenotype2 <- function(genotype)
{
	if (genotype == "AA") return ("B")
	if (genotype == "Aa") return ("O")
	if (genotype == "aa") return ("Y")
	assert("Unknown genotype, should be 'AA', 'Aa' or 'aa'",1==2)
}

GenotypeToPhenotype3 <- function(genotype)
{
	if (genotype == "AA") return ("Y")
	if (genotype == "Aa") return ("B")
	if (genotype == "aa") return ("O")
	assert("Unknown genotype, should be 'AA', 'Aa' or 'aa'",1==2)
}

PhenotypeFrequenciesToGenotypeDensities <- function(phenotypes_frequencies,phenotype_to_genenotype_function)
{
	genotypes_frequencies <- phenotypes_frequencies
	names(genotypes_frequencies) <- phenotype_to_genenotype_function(names(phenotypes_frequencies))
	# Sort the columns
	genotypes_frequencies <- genotypes_frequencies[c("AA","Aa","aa")]
	return (genotypes_frequencies)
}

PhenotypeToGenotype1 <- function(phenotypes)
{
	if (length(phenotypes) == 1)
	{
		phenotype <- phenotypes[1]
		if (phenotype == "O") return ("AA")
		if (phenotype == "Y") return ("Aa")
		if (phenotype == "B") return ("aa")
		assert("Unknown phenotype, should be 'B', 'O' or 'Y'",1==2)
	}
	r <- c()
	for (phenotype in phenotypes)
	{
		r <- c(r,PhenotypeToGenotype1(phenotype))
	}
	return (r)
}

PhenotypeToGenotype2 <- function(phenotypes)
{
	if (length(phenotypes) == 1)
	{
		phenotype <- phenotypes[1]
		if (phenotype == "B") return ("AA")
		if (phenotype == "O") return ("Aa")
		if (phenotype == "Y") return ("aa")
		assert("Unknown phenotype, should be 'B', 'O' or 'Y'",1==2)
	}
	r <- c()
	for (phenotype in phenotypes)
	{
		r <- c(r,PhenotypeToGenotype2(phenotype))
	}
	return (r)
}

PhenotypeToGenotype3 <- function(phenotypes)
{
	if (length(phenotypes) == 1)
	{
		phenotype <- phenotypes[1]
		if (phenotype == "Y") return ("AA")
		if (phenotype == "B") return ("Aa")
		if (phenotype == "O") return ("aa")
		assert("Unknown phenotype, should be 'B', 'O' or 'Y'",1==2)
	}
	r <- c()
	for (phenotype in phenotypes)
	{
		r <- c(r,PhenotypeToGenotype3(phenotype))
	}
	return (r)
}

GetGenotypeToPhenotypeFunction <- function(genotype_to_phenotype_function)
{
  if (genotype_to_phenotype_function == 1) { return(GenotypeToPhenotype1) }
  if (genotype_to_phenotype_function == 2) { return(GenotypeToPhenotype2) }
  if (genotype_to_phenotype_function == 3) { return(GenotypeToPhenotype3) }
	assert("Unknown genotype_to_phenotype_function, should be '1', '2' or '3'",1==2)
}

GetPhenotypeToGenotypeFunction <- function(phenotype_to_genotype_function)
{
  if (phenotype_to_genotype_function == 1) { return(PhenotypeToGenotype1) }
  if (phenotype_to_genotype_function == 2) { return(PhenotypeToGenotype2) }
  if (phenotype_to_genotype_function == 3) { return(PhenotypeToGenotype3) }
	assert("Unknown phenotype_to_genotype_function, should be '1', '2' or '3'",1==2)
}

GetGenotypeToPhenotypeFunctionDescription <- function(genotype_to_phenotype_function)
{
  if (genotype_to_phenotype_function == 1) { return("AA: O, Aa: Y, aa: B") }
  if (genotype_to_phenotype_function == 2) { return("AA: B, Aa: O, aa: Y") }
  if (genotype_to_phenotype_function == 3) { return("AA: Y, Aa: B, aa: O") }
	assert("Unknown genotype_to_phenotype_function, should be '1', '2' or '3'",1==2)
}

PhenotypeToPlotColor <- function(phenotype)
{
	if (phenotype == "Y") return ("yellow")
	if (phenotype == "B") return ("blue")
	if (phenotype == "O") return ("orange")
	assert("Unknown phenotype, should be 'AA', 'Aa' or 'aa'",1==2)
}



RunSimulation <- function(initial_phenotype_densities,function_index,n_generations,year)
{
	genotype_to_phenotype_function <- GetGenotypeToPhenotypeFunction(function_index)
	phenotype_to_genenotype_function <- GetPhenotypeToGenotypeFunction(function_index)
	genotype_to_phenotype_description <- GetGenotypeToPhenotypeFunctionDescription(function_index)

	n_alleles <- 2
	n_genotypes <- 3
	n_phenotypes <- 3
	
	# phenotypes_in_time: keeps track of all frequencies in time
	t_n_rows <- n_generations
	t_n_cols <- 1 + n_alleles + n_genotypes + n_phenotypes
	t <- data.frame(matrix(0,t_n_rows,t_n_cols))
	colnames(t) <- c("t","pA_gametes","pa_gametes","pAA_juveniles","pAa_juveniles","paa_juveniles","pY","pB","pO")
	rownames(t) <- seq(1:n_generations)
	initial_genotypes <- PhenotypeFrequenciesToGenotypeDensities(initial_phenotype_densities,phenotype_to_genenotype_function)
	initial_alleles <- GetAlleleFrequenciesFromGenotypeFrequencies(initial_genotypes)
	t$pA_gametes[1] <- initial_alleles$A
	t$pa_gametes[1] <- initial_alleles$a
	
	# Plot traits in time
	# Start at t=2, because t=1 denotes the initial values
	for (i in c(1:n_generations))
	{
		t$t[i] <- i + 1990 - 1 # -1 due to R arrays have a first index at 1
		
		
		gamete_A <- t$pA_gametes[i]
		gamete_a <- t$pa_gametes[i]
		
		assert("Gamete allele frequenies must sum up to one",abs(gamete_A + gamete_a - 1.0) < 0.0001)
	
		offspring_AA = (gamete_A*gamete_A)
	  offspring_Aa = 2.0 * gamete_A * gamete_a
	  offspring_aa = gamete_a * gamete_a
		offspring_sum_genotypes <- sum(offspring_AA,offspring_Aa,offspring_aa)
	  assert("Offspring genotype frequencies must sum up to one",abs(offspring_sum_genotypes-1.0) < 0.001)
	
		t$pAA_juveniles[i] <- offspring_AA
		t$pAa_juveniles[i] <- offspring_Aa
		t$paa_juveniles[i] <- offspring_aa
	
		offspring_Y <- 0.0
		offspring_B <- 0.0
		offspring_O <- 0.0
		if (genotype_to_phenotype_function("AA") == "Y") {	offspring_Y <- offspring_AA }
		if (genotype_to_phenotype_function("AA") == "B") {	offspring_B <- offspring_AA }
		if (genotype_to_phenotype_function("AA") == "O") {	offspring_O <- offspring_AA }
		if (genotype_to_phenotype_function("Aa") == "Y") {	offspring_Y <- offspring_Aa }
		if (genotype_to_phenotype_function("Aa") == "B") {	offspring_B <- offspring_Aa }
		if (genotype_to_phenotype_function("Aa") == "O") {	offspring_O <- offspring_Aa }
		if (genotype_to_phenotype_function("aa") == "Y") {	offspring_Y <- offspring_aa }
		if (genotype_to_phenotype_function("aa") == "B") {	offspring_B <- offspring_aa }
		if (genotype_to_phenotype_function("aa") == "O") {	offspring_O <- offspring_aa }
		offspring_sum_phenotypes <- sum(offspring_Y,offspring_B,offspring_O)
		assert("Sum of phenotype frequencies must sum to one",
			abs(offspring_sum_phenotypes-1.0) < 0.001
		)
	  t$pY[i] <- offspring_Y
		t$pB[i] <- offspring_B
		t$pO[i] <- offspring_O
	
		offspring_phenotypes <- CreatePhenotypeFrequencies(offspring_Y,offspring_B,offspring_O)
		fitness_matrix <- CreateFitnessMatrix(year)
		
	  phenotype_fitness <- data.matrix(fitness_matrix)  %*% data.matrix(offspring_phenotypes)
	  
		average_fitness  <- (offspring_Y * phenotype_fitness[1,1]) + (offspring_B * phenotype_fitness[2,1]) + (offspring_O * phenotype_fitness[3,1])
	
		# Replicator equations
		adult_Y <- offspring_Y * phenotype_fitness[1,1] / average_fitness
	  adult_B <- offspring_B * phenotype_fitness[2,1] / average_fitness
	  adult_O <- offspring_O * phenotype_fitness[3,1] / average_fitness
	
		assert("Adult phenotypes must sum up to one", sum(adult_Y+adult_B+adult_O-1.0) < 0.0001)
	  adult_genotypes_AA <- 0
	  adult_genotypes_Aa <- 0
	  adult_genotypes_aa <- 0
		if (genotype_to_phenotype_function("AA") == "Y") {	adult_genotypes_AA <- adult_Y }
		if (genotype_to_phenotype_function("AA") == "B") {	adult_genotypes_AA <- adult_B }
		if (genotype_to_phenotype_function("AA") == "O") {	adult_genotypes_AA <- adult_O }
		if (genotype_to_phenotype_function("Aa") == "Y") {	adult_genotypes_Aa <- adult_Y }
		if (genotype_to_phenotype_function("Aa") == "B") {	adult_genotypes_Aa <- adult_B }
		if (genotype_to_phenotype_function("Aa") == "O") {	adult_genotypes_Aa <- adult_O }
		if (genotype_to_phenotype_function("aa") == "Y") {	adult_genotypes_aa <- adult_Y }
		if (genotype_to_phenotype_function("aa") == "B") {	adult_genotypes_aa <- adult_B }
		if (genotype_to_phenotype_function("aa") == "O") {	adult_genotypes_aa <- adult_O }
	
	
	  adult_gametes_A <- adult_genotypes_AA + (0.5 * adult_genotypes_Aa)
	  adult_gametes_a <- adult_genotypes_aa + (0.5 * adult_genotypes_Aa)
	
		assert("Adult phenotypes must sum up to one", sum(adult_gametes_A + adult_gametes_a -1.0) < 0.0001)
	
		if (i < n_generations)
		{
		  t$pA_gametes[i+1] <- adult_gametes_A
		  t$pa_gametes[i+1] <- adult_gametes_a  
		}
	}
  return (t)
}

TestSimulation <- function()
{
  year <- 1996
	function_index <- 1
  initial_phenotype_densities <- CreatePhenotypeFrequencies(0.1,0.6)


	genotype_to_phenotype_function <- GetGenotypeToPhenotypeFunction(function_index)
	phenotype_to_genenotype_function <- GetPhenotypeToGenotypeFunction(function_index)
	genotype_to_phenotype_description <- GetGenotypeToPhenotypeFunctionDescription(function_index)

	initial_genotype_frequencies <- PhenotypeFrequenciesToGenotypeDensities(initial_phenotype_densities,phenotype_to_genenotype_function)
	initial_allele_frequencies <- GetAlleleFrequenciesFromGenotypeFrequencies(initial_genotype_frequencies)
  initial_year <- 1990
		

	# LOOP
	
  current_year <- initial_year
	current_genotype_frequencies <- initial_genotype_frequencies
	current_allele_frequencies <- GetAlleleFrequenciesFromGenotypeFrequencies(current_genotype_frequencies)

		
	t$pAA_juveniles[i] <- offspring_AA
	t$pAa_juveniles[i] <- offspring_Aa
	t$paa_juveniles[i] <- offspring_aa

	offspring_Y <- 0.0
	offspring_B <- 0.0
	offspring_O <- 0.0
	if (genotype_to_phenotype_function("AA") == "Y") {	offspring_Y <- offspring_AA }
	if (genotype_to_phenotype_function("AA") == "B") {	offspring_B <- offspring_AA }
	if (genotype_to_phenotype_function("AA") == "O") {	offspring_O <- offspring_AA }
	if (genotype_to_phenotype_function("Aa") == "Y") {	offspring_Y <- offspring_Aa }
	if (genotype_to_phenotype_function("Aa") == "B") {	offspring_B <- offspring_Aa }
	if (genotype_to_phenotype_function("Aa") == "O") {	offspring_O <- offspring_Aa }
	if (genotype_to_phenotype_function("aa") == "Y") {	offspring_Y <- offspring_aa }
	if (genotype_to_phenotype_function("aa") == "B") {	offspring_B <- offspring_aa }
	if (genotype_to_phenotype_function("aa") == "O") {	offspring_O <- offspring_aa }
	offspring_sum_phenotypes <- sum(offspring_Y,offspring_B,offspring_O)
	assert("Sum of phenotype frequencies must sum to one",
		abs(offspring_sum_phenotypes-1.0) < 0.001
	)
  t$pY[i] <- offspring_Y
	t$pB[i] <- offspring_B
	t$pO[i] <- offspring_O

	offspring_phenotypes <- CreatePhenotypeFrequencies(offspring_Y,offspring_B,offspring_O)
	fitness_matrix <- CreateFitnessMatrix(year)
	
  phenotype_fitness <- data.matrix(fitness_matrix)  %*% data.matrix(offspring_phenotypes)
  
	average_fitness  <- (offspring_Y * phenotype_fitness[1,1]) + (offspring_B * phenotype_fitness[2,1]) + (offspring_O * phenotype_fitness[3,1])

	# Replicator equations
	adult_Y <- offspring_Y * phenotype_fitness[1,1] / average_fitness
  adult_B <- offspring_B * phenotype_fitness[2,1] / average_fitness
  adult_O <- offspring_O * phenotype_fitness[3,1] / average_fitness

	assert("Adult phenotypes must sum up to one", sum(adult_Y+adult_B+adult_O-1.0) < 0.0001)
  adult_genotypes_AA <- 0
  adult_genotypes_Aa <- 0
  adult_genotypes_aa <- 0
	if (genotype_to_phenotype_function("AA") == "Y") {	adult_genotypes_AA <- adult_Y }
	if (genotype_to_phenotype_function("AA") == "B") {	adult_genotypes_AA <- adult_B }
	if (genotype_to_phenotype_function("AA") == "O") {	adult_genotypes_AA <- adult_O }
	if (genotype_to_phenotype_function("Aa") == "Y") {	adult_genotypes_Aa <- adult_Y }
	if (genotype_to_phenotype_function("Aa") == "B") {	adult_genotypes_Aa <- adult_B }
	if (genotype_to_phenotype_function("Aa") == "O") {	adult_genotypes_Aa <- adult_O }
	if (genotype_to_phenotype_function("aa") == "Y") {	adult_genotypes_aa <- adult_Y }
	if (genotype_to_phenotype_function("aa") == "B") {	adult_genotypes_aa <- adult_B }
	if (genotype_to_phenotype_function("aa") == "O") {	adult_genotypes_aa <- adult_O }


  adult_gametes_A <- adult_genotypes_AA + (0.5 * adult_genotypes_Aa)
  adult_gametes_a <- adult_genotypes_aa + (0.5 * adult_genotypes_Aa)

	assert("Adult phenotypes must sum up to one", sum(adult_gametes_A + adult_gametes_a -1.0) < 0.0001)

	if (i < n_generations)
	{
	  t$pA_gametes[i+1] <- adult_gametes_A
	  t$pa_gametes[i+1] <- adult_gametes_a  
	}
	
}