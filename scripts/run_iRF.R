# Predict PPs from TF spatial expression profiles for TFs 
library(iRF)
library(parallel)
library(dplyr)
library(data.table)
library(Matrix)
library(stringr)

# load late data 
n.cores <- 3

# these the PPs for which we will predict 
pp.predict <- c(3, 7)

path <- paste0('./irf_fits/')
dir.create(path, recursive=TRUE)

#################
# load late data
#################
# images
x <- read.csv('../data/clean_x.csv')[, -1]

# dictionary 
dict_file <- './staNMFDicts/K=19/factorization_99.csv'
dict <- read.csv(file = dict_file, header = FALSE)[, -1]

# coefficients
alpha_file <- './staNMFDicts/K=19/alpha_99.csv'
alpha <- read.csv(alpha_file)[, -1]

stopifnot(dim(x)[1] == dim(dict)[1])
stopifnot(dim(x)[2] == dim(alpha)[2])
stopifnot(dim(dict)[2] == dim(alpha)[1])


#####################
# irf parameters
####################
# the threshold to be a gut gene
response <- apply(dict[, pp.predict], 1, max)
thresh.y <- quantile(response, 0.9)

n.iter <- 25
n.bootstrap <- 100
rit.param <- list(depth=5, ntree=2500, 
                  nchild=2, class.id=1, 
                  min.nd=1)

set.seed(47)
n_pixels <- dim(x)[1]
train.id <- sample(1:n_pixels, floor(n_pixels * 0.75))
test.id <- setdiff(1:n_pixels, train.id[[1]])

runReplicate <- function(ii, thresh.y, path, loc, n.cores) {  
  set.seed(ii)

  # keep genes with expression in *any* of these pps
  which_keep_bool <- colSums(alpha[c(3, 4, 6, 7, 11, 15), ]) > 0 # for eye pps
  # which_keep_bool <- colSums(alpha[c(2, 8, 9, 10, 13, 18), ]) > 0 # for gut pps
  gn.keep <- colnames(alpha)[which_keep_bool]
  
  print('genes kept')
  print(length(gn.keep))
  # Convert dictionary atom for given pp into binary response by thresholding.
  y.late <- as.factor(response > thresh.y)
  x.late <- x[,gn.keep]
  
  foo <- tools::toTitleCase(colnames(x.late)) 
  colnames(x.late) <- as.character(1:dim(x.late)[2])

  fit <- iRF(x=x.late[train.id,], 
             y=y.late[train.id], 
             xtest=x.late[test.id,], 
             ytest=y.late[test.id], 
             n.iter=n.iter, 
             n.core=n.cores,
             rit.param=rit.param,
             select.iter=TRUE,
             n.bootstrap=n.bootstrap,
             verbose=TRUE)
  
  filename <- 'irfSpatialFit_foo'
  colnames(x.late) <- foo
  save(file=paste0(path, filename, '.Rdata'), 
       fit, train.id, test.id, x.late, y.late)
}

runReplicate(ii=1, thresh.y=thresh.y, path=path, 
             loc=FALSE, n.cores=n.cores)
