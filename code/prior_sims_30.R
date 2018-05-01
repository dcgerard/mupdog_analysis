## Simulation code for comparing the various priors
## Same as prior_sims.R except only have 50 individuals
library(updog)
library(tidyverse)
library(mupdog)
library(fitPoly)
data("snpdat")
nvec <- rep(100, length = 30)

mean_val <- 4.5

ploidy      <- 6
allele_freq <- mean_val / ploidy
p1geno      <- 4
p2geno      <- 5
bb_od       <- 0.1
mu          <- mean_val
sigma       <- 0.5 ## sd of binomial is about 1.061, so underdispersed here.

pilist <- list()
pilist[[1]] <- dbinom(x = 0:ploidy, size = ploidy, prob = allele_freq)
pilist[[2]] <- dbetabinom(x = 0:ploidy, size = ploidy, mu = allele_freq, rho = bb_od, log = FALSE)
pilist[[3]] <- dnorm(x = 0:ploidy, mean = mu, sd = sigma)
pilist[[3]] <- pilist[[3]] / sum(pilist[[3]])
pilist[[4]] <- 0:ploidy
pilist[[4]] <- pilist[[4]] / sum(pilist[[4]])
pilist[[5]] <- updog::get_q_array(ploidy = ploidy)[p1geno + 1, p2geno + 1, ]
pilist[[6]] <- c(1, 4, 1, 4, 1, 4, 8)
pilist[[6]] <- pilist[[6]] / sum(pilist[[6]])
pilist[[7]] <- rep(1 / (ploidy + 1), length = ploidy + 1)
pinames <- c("hw", "bb", "norm", "ash", "f1", "flex", "uniform")
names(pilist) <- pinames

saveRDS(object = pilist, file = "../output/prior_sims_30/pilist.RDS")

seq <- 0.001
odvec <- c(0, 0.005, 0.01)
biasvec <- c(1, 0.75, 0.5)
itermax <- 500

parvals <- expand.grid(seed = 1:itermax, od = odvec, bias = biasvec, geno_dist = pinames)
parvals$seq <- seq
parvals$ploidy <- ploidy
set.seed(1)
parvals <- parvals[sample(1:nrow(parvals)), ] ## randomize order so heavy computation doesn't cluster together
parlist <- split(parvals, seq(nrow(parvals)))

one_rep <- function(args, nvec, pilist) {
  set.seed(args$seed)
  pivec <- pilist[[args$geno_dist]]
  geno <- mupdog::rgeno(n = length(nvec), ploidy = args$ploidy, model = "flex", pivec = pivec)
  refvec <- mupdog::rflexdog(sizevec = nvec, geno = geno, ploidy = args$ploidy, seq = args$seq, bias = args$bias, od = args$od)

  bias_init_vec <- exp(c(-0.7, -0.3, 0, 0.3, 0.7))
  lbest_vec <- rep(-Inf, length = 7)
  mout <- list()
  mout_temp <-  list()
  for (fit_index in 1:length(bias_init_vec)) {
    ## Fit methods
    mout_temp[[1]] <- mupdog::flexdog(refvec = refvec, sizevec = nvec, ploidy = args$ploidy, model = "hw", verbose = FALSE, bias = bias_init_vec[fit_index])
    mout_temp[[2]] <- mupdog::flexdog(refvec = refvec, sizevec = nvec, ploidy = args$ploidy, model = "bb", verbose = FALSE, bias = bias_init_vec[fit_index])
    mout_temp[[3]] <- mupdog::flexdog(refvec = refvec, sizevec = nvec, ploidy = args$ploidy, model = "norm", verbose = FALSE, bias = bias_init_vec[fit_index])
    mout_temp[[4]] <- mupdog::flexdog(refvec = refvec, sizevec = nvec, ploidy = args$ploidy, model = "ash", verbose = FALSE, bias = bias_init_vec[fit_index])
    mout_temp[[5]] <- mupdog::flexdog(refvec = refvec, sizevec = nvec, ploidy = args$ploidy, model = "f1", verbose = FALSE, bias = bias_init_vec[fit_index])
    mout_temp[[6]] <- mupdog::flexdog(refvec = refvec, sizevec = nvec, ploidy = args$ploidy, model = "flex", verbose = FALSE, bias = bias_init_vec[fit_index])
    suppressWarnings(mout_temp[[7]] <- mupdog::flexdog(refvec = refvec, sizevec = nvec, ploidy = args$ploidy, model = "uniform", verbose = FALSE, bias = bias_init_vec[fit_index]))

    ## Choose highest likelihood
    for (mindex in 1:length(mout_temp)) {
     if (mout_temp[[mindex]]$llike > lbest_vec[mindex]) {
       mout[[mindex]] <- mout_temp[[mindex]]
       lbest_vec[mindex] <- mout_temp[[mindex]]$l
     }
    }
  }

  ## Compare to fitPoly
  ## Won't run with too few samples
  # tempdf <- data.frame(MarkerName = "SNP",
  #                      SampleName = 1:length(nvec),
  #                      ratio = refvec / nvec)
  # fpout <- fitPoly::fitOneMarker(ploidy = ploidy,
  #                                marker = "SNP",
  #                                data = tempdf)

  ## Save output
  pc_vec <- rep(NA, length = length(mout))
  names(pc_vec) <- paste0("pc_", names(pilist))

  epm_vec <- rep(NA, length = length(mout))
  names(epm_vec) <- paste0("epm_", names(pilist))

  seq_vec <- rep(NA, length = length(mout))
  names(seq_vec) <- paste0("seq_", names(pilist))

  bias_vec <- rep(NA, length = length(mout))
  names(bias_vec) <- paste0("bias_", names(pilist))

  od_vec <- rep(NA, length = length(mout))
  names(od_vec) <- paste0("od_", names(pilist))

  for (index in 1:length(mout)) {
    pc_vec[index]   <- mean(mout[[index]]$geno == geno)
    seq_vec[index]  <- mout[[index]]$seq
    bias_vec[index] <- mout[[index]]$bias
    od_vec[index]   <- mout[[index]]$od
    epm_vec[index]  <- mout[[index]]$prop_mis
  }

  return_vec <- c(pc_vec, epm_vec, seq_vec, bias_vec, od_vec)

  return(return_vec)
}


library(parallel)
nc <- 6
cl   <- makeCluster(nc)
cat("Running multithreaded computations with", nc, "threads.\n")
sout <- t(parallel::parSapply(cl = cl, X = parlist, FUN = one_rep,
                              nvec = nvec, pilist = pilist))
stopCluster(cl)

outtab <- cbind(parvals, sout)
saveRDS(object = outtab, file = "../output/prior_sims_30/sims_out.RDS")


