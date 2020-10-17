 # Write your own function for calculating the sums of squares of the whole dataset. SS(x) is the sum of the squared deviations from the mean given by:  
#  $\sum (x_i- \bar{x})^2$  

my_ss <- function(v) {sum((v - mean(v))^2)}


# Put the code for the function in its own script and call it from your main script
file <- here::here("functions", "my_ss.R")
source(file)
my_ss(chaff2$mass)
