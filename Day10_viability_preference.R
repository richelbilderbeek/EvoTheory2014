rm(list = ls())
options(error = browser)
options(echo = FALSE)

Vz <- 0.2
Vy <- Vz
Cxy <- 0.1
theta_male <- 10 # Male optimum value for viability selection
omega_male <- 1  #Distribution of viability value around theta_male
theta_female <- 5 # Female optimum value for viability selection
omega_female <- 1  #Distribution of viability value around theta_female

n_generations <- 100
initial_mean_z <- 0 #Male trait
initial_mean_y <- 0 #Female preference

male_color <- "red"
female_color <- "blue"

phenotypes_n_rows <- 2
phenotypes_n_cols <- 1
phenotypes <- data.frame(matrix(0,phenotypes_n_rows,phenotypes_n_cols))
phenotypes[1,1] <- "male trait"
phenotypes[2,1] <- "female preference"
rownames(phenotypes) <- c("z","y")
colnames(phenotypes) <- c("average")

# phenotypes_in_time: keeps track of all frequencies in time
t_n_rows <- n_generations
t_n_cols <- phenotypes_n_rows
phenotypes_in_time <- data.frame(matrix(0,t_n_rows,t_n_cols))
colnames(phenotypes_in_time) <- rownames(phenotypes)
rownames(phenotypes_in_time) <- seq(1:n_generations)
phenotypes_in_time[1,1] <- initial_mean_z
phenotypes_in_time[1,2] <- initial_mean_y

# Plot traits in time
# Start at t=2, because t=1 denotes the initial values
for (t in c(2:n_generations))
{
	cur_z <- phenotypes_in_time[t-1,1]
	cur_y <- phenotypes_in_time[t-1,2]
  beta_z <- cur_y - (cur_z - theta_male) / (omega_male * omega_male)
  beta_y <- -(cur_y - theta_female) / (omega_female * omega_female)
	delta_z <- 0.5 * ((Vz  * beta_z) + (Vz  * beta_y))
	delta_y <- 0.5 * ((Cxy * beta_z) + (Vy  * beta_y))
	
	next_z <- cur_z + delta_z
	next_y <- cur_y + delta_y
	
	phenotypes_in_time[t,1] <- next_z
	phenotypes_in_time[t,2] <- next_y
}

# Need to use as.matrix, otherwise 'graphical parameter "type" is obsolete'
plot(
	as.matrix(phenotypes_in_time[1]),
	xlim=c(1,n_generations),
	ylim=c(min(phenotypes_in_time),max(phenotypes_in_time)),
	type="l",col=male_color,
	main=paste("Average phenotypes in time\nfor Vz=",Vz,", Cxy=",Cxy,",\ntheta_male=",theta_male,", omega_male=",omega_male,",\ntheta_female=",theta_female,", omega_female=",omega_female),
	ylab="Trait value",
	xlab="Time (generations)"
)
legend_x <- 0.6 * n_generations
legend_y <- max(phenotypes_in_time) - (0.4 * (max(phenotypes_in_time) - min(phenotypes_in_time)))
legend(
	legend_x,
	legend_y,
	c("Male trait","Female preference"),col=c(male_color,female_color),pch=c(16,16)
)
lines(phenotypes_in_time[2],col=female_color)
lines(c(1,n_generations),c(theta_male,theta_male),lty="dashed",col=male_color)
lines(c(1,n_generations),c(theta_female,theta_female),lty="dashed",col=female_color)


# Phase plot
min_z <- -20
max_z <-  20
res_z <- 20 #Resolution: number of arrows in horizontal direction
min_y <- -20
max_y <-  20
res_y <- 20 #Resolution: number of arrows in vertical direction
plot(
	"",
	xlim=c(min_z,max_z),
	ylim=c(min_y,max_y),
	type="n",
	main=paste("Phenotype change\nfor Vz=",Vz,", Cxy=",Cxy,",\ntheta_male=",theta_male,", omega_male=",omega_male,",\ntheta_female=",theta_female,", omega_female=",omega_female),
	ylab="Female preference (y)",
	xlab="Male trait value (z)"
)
equilibrium_line <- function(z) 
{
  y <- (z - theta_male) / (omega_male * omega_male)	
}
lines(c(min_z,max_z),c(equilibrium_line(min_z),equilibrium_line(max_z)),lty="dashed",col="black")
lines(c(theta_male,theta_male),c(min_z,max_z),lty="dashed",col=male_color)
lines(c(min_z,max_z),c(theta_female,theta_female),lty="dashed",col=female_color)

for (cur_z in seq(min_z,max_z,(max_z-min_z)/res_z))
{
	for (cur_y in seq(min_y,max_y,(max_y-min_y)/res_y))
	{
	  beta_z <- cur_y - (cur_z - theta_male) / (omega_male * omega_male)
	  beta_y <- -(cur_y - theta_female) / (omega_female * omega_female)
		delta_z <- 0.5 * ((Vz  * beta_z) + (Vz  * beta_y))
		delta_y <- 0.5 * ((Cxy * beta_z) + (Vy  * beta_y))
    arrows(
    	cur_z,
    	cur_y,
    	cur_z + (0.5 * delta_z),
    	cur_y + (0.5 * delta_y),
    	length=0.05
    )
	}
}