source('../utils/olsLassoFit.R')
load('../data/lateData.Rdata')

x <- read.csv('../data/clean_x.csv')[, -1]

results_dir <- './staNMFDicts/'

args <- commandArgs(trailingOnly = TRUE)

K <- as.integer(args[1])

print(paste0('getting PP coeffs for fit with K = ', K))

get_PP_coeffs <- function(x, dict){
  # takes in images in the columns of x and fits the PP coefficients
  # using lasso to dict
  PP_coeffs <- matrix(0, dim(dict)[2], dim(x)[2])
  for(i in 1:dim(x)[2]){
    PP_coeffs[, i] <- fitBeta(x[, i], dict)
    
    if((i %% 20) == 0){
      prnt <- paste0('(', as.character(i), '/', as.character(dim(x)[2]),
                     ') images done')
      print(prnt)
    }
  }
  colnames(PP_coeffs) <- colnames(x)
  return(PP_coeffs)
}

# load dictionary 
# for now ... we pick any restart
dict_file_k <- paste0(results_dir, 'K=', K, '/factorization_99.csv')
dict <- read.csv(file = dict_file_k, header = FALSE)[, -1]

alpha <- get_PP_coeffs(x, as.matrix(dict))

outfile <- paste0(results_dir, 'K=', K, '/alpha_99.csv')
print(paste0('saving fit into: ', outfile))
write.csv(alpha, file = outfile)
