rm(list = ls())
options(error = browser)
options(echo = FALSE)
library(klaR)
library(testit)

CreateFitnessMatrix <- function()
{
	fitness_matrix_n_rows <- 3
	fitness_matrix_n_cols <- 3
	fitness_matrix <- data.frame(matrix(0,fitness_matrix_n_rows,fitness_matrix_n_cols))
	fitness_matrix[1,1] <- 1.0
	fitness_matrix[2,1] <- 2.0
	fitness_matrix[3,1] <- 0.3
	fitness_matrix[1,2] <- 1.0
	fitness_matrix[2,2] <- 1.0
	fitness_matrix[3,2] <- 4.0
	fitness_matrix[1,3] <- 2.0
	fitness_matrix[2,3] <- 0.1
	fitness_matrix[3,3] <- 1.0
	rownames(fitness_matrix) <- c("Y","B","O")
	colnames(fitness_matrix) <- c("Y","B","O")
	return (fitness_matrix)
}

CreatePhenotypeFrequencies <- function(pY = 0.0, pB = 0.0, pO = 0.0)
{
	n_rows <- 3
	n_cols <- 1
	fs <- data.frame(matrix(0,n_rows,n_cols))
	fs[1,1] <- pY
	fs[2,1] <- pB
	fs[3,1] <- pO
	colnames(fs) <- c("f")
	rownames(fs) <- c("Y","B","O")
	assert("Phenotype Frequencies must sum up to one or be empty",
		abs(sum(fs) - 1.0) < 0.0001
		|| (pY == 0.0 && pB == 0.0 && pO == 0.0)
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

PhenotypeToPlotColor <- function(phenotype)
{
	if (phenotype == "Y") return ("yellow")
	if (phenotype == "B") return ("blue")
	if (phenotype == "O") return ("orange")
	assert("Unknown phenotype, should be 'AA', 'Aa' or 'aa'",1==2)
}



RunSimulation <- function(initial_A,genotype_to_phenotype_function,n_generations)
{
	n_alleles <- 2
	n_genotypes <- 3
	n_phenotypes <- 3
	
	# phenotypes_in_time: keeps track of all frequencies in time
	t_n_rows <- n_generations
	t_n_cols <- 1 + n_alleles + n_genotypes + n_phenotypes
	t <- data.frame(matrix(0,t_n_rows,t_n_cols))
	colnames(t) <- c("A","a","AA","Aa","aa","Y","B","O")
	rownames(t) <- seq(1:n_generations)
	t[1,1] = 1
	t[1,2] = initial_A
	t[1,3] = 1.0 - initial_A
	t
	
	# Plot traits in time
	# Start at t=2, because t=1 denotes the initial values
	for (i in c(2:n_generations))
	{
		t[i,1] <- i
		gamete_A <- t[i-1,2]
		gamete_a <- t[i-1,3]
		assert("Gamete allele frequenies must sum up to one",abs(gamete_A + gamete_a - 1.0) < 0.0001)
	
		offspring_AA = (gamete_A*gamete_A)
	  offspring_Aa = 2.0 * gamete_A * gamete_a
	  offspring_aa = gamete_a * gamete_a
		offspring_sum_genotypes <- sum(offspring_AA,offspring_Aa,offspring_aa)
	  assert("Offspring genotype frequencies must sum up to one",abs(offspring_sum_genotypes-1.0) < 0.001)
	
		t[i-1,4] <- offspring_AA
		t[i-1,5] <- offspring_Aa
		t[i-1,6] <- offspring_aa
	
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
	  t[i-1,7] <- offspring_Y
		t[i-1,8] <- offspring_B
		t[i-1,9] <- offspring_O
	
		offspring_phenotypes <- CreatePhenotypeFrequencies(offspring_Y,offspring_B,offspring_O)
		fitness_matrix <- CreateFitnessMatrix()
	  assert("Must work for matrix multiplication",nrow(fitness_matrix) == ncol(fitness_matrix))
	
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
	
		t[i,2] <- adult_gametes_A
		t[i,3] <- adult_gametes_a  
	}
	t <- t[-nrow(t),] 
	
	plot(
		as.matrix(t[2]),ylim=c(0,1),
		type="l",col="red",
		main="Allele frequencies in time",
		xlab="Time (generations)",ylab="Allele frequency"
	)
	lines(as.matrix(t[3]),col="blue")
	legend_x <- 0.6 * n_generations
	legend_y <- 1.0
	legend(
		legend_x,
		legend_y,
		c("A","a"),col=c("red","blue"),pch=c(16,16)
	)
	
	plot(
		as.matrix(t[4]),ylim=c(0,1),
		type="l",col=PhenotypeToPlotColor(genotype_to_phenotype_function("AA")),
		main="Genotype frequencies in time",
		xlab="Time (generations)",ylab="Genotype frequency"
	)
	lines(as.matrix(t[5]),col=PhenotypeToPlotColor(genotype_to_phenotype_function("Aa")))
	lines(as.matrix(t[6]),col=PhenotypeToPlotColor(genotype_to_phenotype_function("aa")))
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

	triplot(as.matrix(t[,c(4:6)]),
	  label=c(genotype_to_phenotype_function("AA"),genotype_to_phenotype_function("Aa"),genotype_to_phenotype_function("aa")),
	  grid = TRUE,
		type="l"
	)

  return (t)
}


genotype_to_phenotype_function <- GenotypeToPhenotype2
initial_A <- 0.5
n_generations <- 50
t <- RunSimulation(initial_A,genotype_to_phenotype_function,n_generations)

