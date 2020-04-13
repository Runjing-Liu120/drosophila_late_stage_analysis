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
pp.predict <- c(8)

path <- paste0('./irf_fits_2020-04-13/')
dir.create(path, recursive=TRUE)
filename <- 'irfSpatialFit_gut_set12pp'

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
  p.sample <- 1
  set.seed(53453)
  filename <- filename
}else{
  p.sample <- as.numeric(args[1])
  set.seed(as.numeric(args[2]))
  filename <- args[3]
}

#################
# load late data
#################
# take care of replicates
CombineReplicates <- function(x){
  gene.names.clean <- gsub('\\..', '', colnames(x))
  unique.genes <- unique(gene.names.clean)
  length(unique.genes)
  x.clean <- matrix(0, nrow = dim(x)[1], ncol = length(unique.genes))
  colnames(x.clean) <- unique.genes
  for(i in 1:length(unique.genes)){
    gene.i <- unique.genes[i]
    x.gene <- x[, gene.names.clean == gene.i]
    if(is.null(dim(x.gene))){
      stopifnot(sum(gene.names.clean == gene.i) == 1)
      x.clean[, gene.i] <- x.gene
    }else{
      stopifnot(sum(gene.names.clean == gene.i) > 1)
      x.clean[, gene.i] <- rowMeans(x.gene)
    }
  }
  return(x.clean)
}

# images
x.raw <- read.csv('../data/clean_x.csv')[, -1]
x <- CombineReplicates(x.raw)
# x <- x[, colnames(x) != 'Pdp1']

# dictionary 
dict_file <- './staNMFDicts/K=12/factorization_99.csv'
dict <- read.csv(file = dict_file, header = FALSE)[, -1]

# coefficients
alpha_file <- './staNMFDicts/K=12/alpha_99.csv'
alpha.raw <- read.csv(alpha_file)[, -1]
alpha <- CombineReplicates(alpha.raw)
# alpha <- alpha[, colnames(alpha) != 'Pdp1']

# subsample 
if(p.sample < 1){
  n.genes <- dim(x)[2]
  n.sample <- round(n.genes * p.sample)
  genes.sampled <- sample(n.genes, n.sample)
  x <- x[, genes.sampled]
  alpha <- alpha[, genes.sampled]
}

stopifnot(dim(x)[1] == dim(dict)[1])
stopifnot(dim(x)[2] == dim(alpha)[2])
stopifnot(dim(dict)[2] == dim(alpha)[1])

#####################
# irf parameters
####################
# the threshold to be a gut gene
# response <- apply(dict[, pp.predict], 1, max)
response <- dict[, pp.predict]
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
  # which_keep_bool <- colSums(alpha[c(3, 4, 6, 7, 11, 15), ]) > 0 # for eye pps
  # which_keep_bool <- colSums(alpha[c(2, 8, 9, 10, 13, 18), ]) > 0 # for gut pps
  # which_keep_bool <- colSums(alpha[c(10, 11), ]) > 0 # for eye pps
  which_keep_bool <- colSums(alpha[c(4, 6, 8, 12), ]) > 0 # for gut pps
  gn.keep <- colnames(alpha)[which_keep_bool]
  
  print('genes kept')
  print(length(gn.keep))
  # Convert dictionary atom for given pp into binary response by thresholding.
  y.late <- as.factor(response > thresh.y)
  x.late <- x[,gn.keep]
  
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
  
  outfile <- paste0(path, filename, '.Rdata')
  print('saving iRF results into: ')
  print(outfile)
  save(file=outfile, 
       fit, train.id, test.id, x.late, y.late)
}

runReplicate(ii=1, thresh.y=thresh.y, path=path, 
             loc=FALSE, n.cores=n.cores)
