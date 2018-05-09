library(iterators)
library(updog)
library(parallel)
library(foreach)
library(doParallel)
ploidy_seq <- c(2, 4, 6)
seq        <- 0.001
bias_seq   <- seq(0.5, 1, length = 11)
od_seq     <- seq(0, 0.02, length = 11)
nseq       <- seq_len(1000)
alpha_seq  <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
as.data.frame(expand.grid(n      = nseq,
                          ploidy = ploidy_seq,
                          seq    = seq,
                          bias   = bias_seq,
                          od     = od_seq,
                          alpha  = alpha_seq)) ->
  odat

## so that computation is not heavy in some clusters
set.seed(1)
odat_order <- sample(seq_len(nrow(odat)))
odat <- odat[odat_order, ]


nc <- 6
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl = cl)
stopifnot(foreach::getDoParWorkers() > 1)

glist <- foreach(index = seq_len(nrow(odat)),
                 .export = c("oracle_joint")) %dopar% {
  dist <- stats::dbinom(x = 0:odat$ploidy[index],
                        size = odat$ploidy[index],
                        prob = odat$alpha[index])
  oracle_joint(n = odat$n[index],
               ploidy = odat$ploidy[index],
               seq = odat$seq[index],
               bias = odat$bias[index],
               od = odat$od[index],
               dist = dist)
}
stopCluster(cl)
odat$joint <- glist

odat$pmiss <- sapply(X = odat$joint, FUN = oracle_mis_from_joint)
odat$cor   <- sapply(X = odat$joint, FUN = oracle_cor_from_joint)
saveRDS(object = odat, file = "../output/oracle_explore/odat.RDS")
