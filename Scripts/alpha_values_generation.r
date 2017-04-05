setwd('~/MSP_Model/Output/Data/')
install.packages('gtools', repos = 'https://cran.mtu.edu/')
require(gtools)
epsilon <- .20 # Epsilon step size, default is 0.20
a_values <- seq(from = 0, to = 1, by = epsilon) # The unique values for each sector and site
alpha <- permutations(n = length(a_values),7,a_values,repeats.allowed=T) # Matrix of unique alpha weights
write.csv(alpha, file = file.path(paste0(getwd(),'/alpha_weights.csv')))
system(paste('open ',file.path(paste0(getwd(),'/alpha_weights.csv'))))
