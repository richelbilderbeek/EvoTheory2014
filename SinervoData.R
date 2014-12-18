library(testit)
library(klaR)

CreateFitnessMatrix1996 <- function()
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

CreateFitnessMatrix2001 <- function()
{
	fitness_matrix_n_rows <- 3
	fitness_matrix_n_cols <- 3
	fitness_matrix <- data.frame(matrix(0,fitness_matrix_n_rows,fitness_matrix_n_cols))
	fitness_matrix[1,1] <- 1.0
	fitness_matrix[2,1] <- 2.0
	fitness_matrix[3,1] <- 1.2
	fitness_matrix[1,2] <- 0.7
	fitness_matrix[2,2] <- 1.0
	fitness_matrix[3,2] <- 1.8
	fitness_matrix[1,3] <- 2.2
	fitness_matrix[2,3] <- 0.9
	fitness_matrix[3,3] <- 1.0
	rownames(fitness_matrix) <- c("Y","B","O")
	colnames(fitness_matrix) <- c("Y","B","O")
	return (fitness_matrix)
}

CreateFitnessMatrix <- function(year)
{
  if (year == 1996) return(CreateFitnessMatrix1996())	
  if (year == 2001) return(CreateFitnessMatrix2001())	
	assert("Unknown year, should be '1996' or '2001'",1==2)
}

CreatePredictedPhenotypesRaw <- function()
{
	t <- data.frame(
		matrix(
		  c( 
				0.33,0.57,0.11,
				0.13,0.75,0.12,
				0.36,0.33,0.31,
				0.42,0.32,0.26,
				0.49,0.39,0.13,
				0.35,0.55,0.10,
				0.43,0.45,0.12,
				0.60,0.15,0.25,
				0.52,0.37,0.10,
				0.45,0.46,0.09
		  ),
			10,
			3,
		  byrow = TRUE
		)
	)
	names(t) = c("Y","B","O")
	rownames(t) = seq(1990,1999)
  return (t)
}	

CreatePredictedPhenotypes <- function()
{
  t <- CreatePredictedPhenotypesRaw()
	return (t / rowSums(t))
}

PlotPredictedPhenotypes <- function()
{
	rowSums(CreatePredictedPhenotypes())
	triplot(as.matrix(CreatePredictedPhenotypes()),
		main=paste("Observed phenotypes in time\n(Sinervo & Lively (1996), Nature)"),
	  label=c("Y","B","O"),
	  grid = TRUE,
		type="l",
		col="black"
	)
}

PlotPredictedPhenotypesInTime <- function()
{
	png("SinervoPhenotypesInTime.png")
	plot(
		CreatePredictedPhenotypes()$Y,
		main=paste("Observed phenotypes in time\n(Sinervo & Lively (1996), Nature)"),
		type="l",
		ylim=c(0,1),
		ylab="Frequency",
		col="yellow",
		lwd = 3
	)
	lines(
		CreatePredictedPhenotypes()$B,
		col="blue",
		lwd = 1
	)
	lines(
		CreatePredictedPhenotypes()$O,
		col="orange",
		lwd = 2
	)
	legend_x <- 4
	legend_y <- 1.0
	legend(
		legend_x,
		legend_y,
		c("yellow","blue","orange"),
		col=c("yellow","blue","orange"),
		pch=c(16,16),
		lwd = c(3,1,2)
	)
	dev.off()
}

PlotPredictedPhenotypesInTime()
