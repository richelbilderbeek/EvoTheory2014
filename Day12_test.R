library(testit)

TestDay12 <- function()
{
	for (function_index in seq(1:3))
	{
		genotype_to_phenotype_function <- GetGenotypeToPhenotypeFunction(function_index)
		phenotype_to_genotype_function <- GetPhenotypeToGenotypeFunction(function_index)
		for (genotype in c("AA","Aa","aa"))
		{
			phenotype <- genotype_to_phenotype_function(genotype)
			genotype_again <- phenotype_to_genotype_function(phenotype)
			assert("",genotype_again == genotype)			
		}
		for (phenotype in c("Y","B","O"))
		{
			genotype <- phenotype_to_genotype_function(phenotype)
			phenotype_again <- genotype_to_phenotype_function(genotype)
			assert("",phenotype_again == phenotype)			
		}
	}
	# PhenotypeToGenotype1: column handling
  {
    assert("",length(PhenotypeToGenotype1(c("O","O"))) == 2)
    assert("",PhenotypeToGenotype1(c("O","O")) == c("AA","AA"))
    assert("",PhenotypeToGenotype1(c("O","O")) != c("Aa","Aa"))
    assert("",PhenotypeToGenotype1(c("O","Y","B")) == c("AA","Aa","aa"))
  }
  # PhenotypeFrequenciesToGenotypeDensities
	{
		phenotype_to_genotype_function <- GetPhenotypeToGenotypeFunction(1)
	  initial_phenotypes <- CreatePredictedPhenotypes()["1990",]
		# 		 Aa        aa        AA
		# 1990 0.3267327 0.5643564 0.1089109

		# For function 1:
		# if (phenotype == "O") return ("AA")
		# if (phenotype == "Y") return ("Aa")
		# if (phenotype == "B") return ("aa")
		initial_genotypes <- PhenotypeFrequenciesToGenotypeDensities(initial_phenotypes,GetPhenotypeToGenotypeFunction(1))
    assert("",initial_phenotypes$O == initial_genotypes$AA)
    assert("",initial_phenotypes$Y == initial_genotypes$Aa)
    assert("",initial_phenotypes$B == initial_genotypes$aa)

		# For function 2:
		# if (phenotype == "B") return ("AA")
		# if (phenotype == "O") return ("Aa")
		# if (phenotype == "Y") return ("aa")
		initial_genotypes <- PhenotypeFrequenciesToGenotypeDensities(initial_phenotypes,GetPhenotypeToGenotypeFunction(2))
    assert("",initial_phenotypes$B == initial_genotypes$AA)
    assert("",initial_phenotypes$O == initial_genotypes$Aa)
    assert("",initial_phenotypes$Y == initial_genotypes$aa)

		# For function 3:
		# if (phenotype == "Y") return ("AA")
		# if (phenotype == "B") return ("Aa")
		# if (phenotype == "O") return ("aa")
		initial_genotypes <- PhenotypeFrequenciesToGenotypeDensities(initial_phenotypes,GetPhenotypeToGenotypeFunction(3))
    assert("",initial_phenotypes$Y == initial_genotypes$AA)
    assert("",initial_phenotypes$B == initial_genotypes$Aa)
    assert("",initial_phenotypes$O == initial_genotypes$aa)
		
	}
}

TestDay12()