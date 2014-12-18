library(testit)

source('SinervoData.R')
source('../R/Genetics.R')

# Creates all fitness matrices as one big table
#
# 1 1 1
# 1 1 1
# 1 1 1
# 2 2 2
# 2 2 2
# 2 2 2
# 3 3 3
# . . .
# CreateAllFitnessMatrices <- function()
# {
# 
# 	 fitness_matrices <- c(
# 		CreateFitnessMatrix1(1996),
# 		CreateFitnessMatrix2(1996),
# 		CreateFitnessMatrix3(1996),
# 		CreateFitnessMatrix1(2001),
# 		CreateFitnessMatrix2(2001),
# 		CreateFitnessMatrix3(2001)
# 	)
# 	for (i in c(1:6))
# 	{
# 		fitness_matrix <- fitness_matrices[i]
# 		assert("",ncol(fitness_matrix) == nrow(fitness_matrix))
# 	}
# 	return (c)
# }

CreateFitnessMatrix1 <- function(year)
{
  m <- CreateFitnessMatrix(year)
	assert("",ncol(m) > 1)
	assert("",nrow(m) > 1)
	colnames(m) <- c("Aa","aa","AA")
	rownames(m) <- c("Aa","aa","AA")

	m <- SortFitnessMatrix(m)

	return (m)	
}

CreateFitnessMatrix2 <- function(year)
{
  m <- CreateFitnessMatrix(year)
	assert("",ncol(m) > 1)
	assert("",nrow(m) > 1)
	colnames(m) <- c("aa","AA","Aa")
	rownames(m) <- c("aa","AA","Aa")

	m <- SortFitnessMatrix(m)

	return (m)	
}

CreateFitnessMatrix3 <- function(year)
{
  m <- CreateFitnessMatrix(year)
	assert("",ncol(m) > 1)
	assert("",nrow(m) > 1)
  colnames(m) <- c("AA","Aa","aa")
	rownames(m) <- c("AA","Aa","aa")

	m <- SortFitnessMatrix(m)

	return (m)	
}

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

RunSimulation <- function(initial_genotype_frequencies_adults,fitness_matrix,n_generations)
{
# 	print("RunSimulation:Test parameters")
# 	initial_genotype_frequencies_adults <- CreateGenotypeFrequencies(0.1,0.2,0.7)
# 	fitness_matrix <- CreateFitnessMatrix1(1996)
# 	n_generations <- 5
# 	initial_genotype_frequencies_adults = CreateGenotypeFrequencies(0.2,0.3,0.5)
# 	print("~RunSimulation:Test parameters")

	assert("",ncol(initial_genotype_frequencies_adults)==ncol(fitness_matrix))
  SortGenotypeFrequencies(initial_genotype_frequencies_adults)
	
	t <- data.frame(matrix(0,0,length(initial_genotype_frequencies_adults)))
	colnames(t) <- colnames(initial_genotype_frequencies_adults)
	
	t <- rbind(initial_genotype_frequencies_adults)

	assert("",as.matrix(t[1,]) == as.matrix(initial_genotype_frequencies_adults))
  for (i in c(2:n_generations))
  {
  	genotype_densities_adults <- t[i-1,]
  	# genotype_densities_adults
    genotype_densities_new_adults <- CalculateNextGenotypeDensities(genotype_densities_adults,fitness_matrix)
  	# genotype_densities_new_adults
  	t <- rbind(t,genotype_densities_new_adults)
  }
	return (t)
}

SortFitnessMatrix <-function(m)
{
	m
	m <- m[c("AA","Aa","aa")]
	m <- m[order(names(m)),]
  
	#assert("",names(m) == "AA","Aa","aa")
	#assert("",rownames(m) == "AA","Aa","aa")
  return (m)
}

SortGenotypes <-function(m)
{
	m <- m[c("AA","Aa","aa")]
  return (m)
}


CalculateNextGenotypeDensities <- function(genotype_frequencies_parents,fitness_matrix)
{
# 	print("CalculateNextGenotypeDensities: Testing arguments")
#   genotype_frequencies_parents <- CreateGenotypeFrequencies(0.2,0.3,0.5)
#   fitness_matrix <- CreateFitnessMatrix1(1996)
# 	print("~CalculateNextGenotypeDensities: Testing arguments")

	assert("",nrow(genotype_frequencies_parents)==1)
	assert("",ncol(genotype_frequencies_parents)==ncol(fitness_matrix))
	
	# Let the genotypes create gametes
	allele_frequencies_gametes <- GetAlleleFrequenciesFromGenotypeFrequencies(genotype_frequencies_parents)

	# Mix the gametes and create juveniles
	genotype_frequencies_juveniles <- GetGenotypeFrequenciesFromAlleleFrequencies(allele_frequencies_gametes)

	# Calculate the genotype fitnesses	
  genotype_fitness_adults <- t(data.matrix(fitness_matrix)  %*% t(data.matrix(genotype_frequencies_juveniles)))
  
	# Normalized by w_bar
	mean_weightes_fitness <- sum(genotype_frequencies_juveniles * genotype_fitness_adults)

	# Replicator equation
	genotype_frequencies_adults <- genotype_frequencies_juveniles * genotype_fitness_adults / mean_weightes_fitness

	assert("",nrow(genotype_frequencies_adults)==1)
  assert("",abs(sum(genotype_frequencies_adults) - 1) < 0.001)
	
  return 	(genotype_frequencies_adults)
}