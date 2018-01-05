## Select toy dataset to include in mupdog package
library(mupdog)
set.seed(1)
n <- 10
p <- 100

mout <- readRDS(file = "../output/uit_fit/fit_uit1.RDS")

good_bias <- mout$bias > 0.5 & mout$bias < 2
good_seq  <- mout$seq < 0.05
good_od   <- mout$od < 0.05
good_allele <- mout$allele_freq < 0.9 & mout$allele_freq > 0.1

good_snp <- good_bias & good_seq & good_od & good_allele


good_ind <- mout$inbreeding < 0.3

which_snp <- sample((1:length(good_snp))[good_snp], p)
which_ind <- sample((1:length(good_ind))[good_ind], n)

refmat  <- mout$input$refmat[which_ind, which_snp]
sizemat <- mout$input$sizemat[which_ind, which_snp]
ploidy  <- mout$input$ploidy

uitdewilligen <- list()
uitdewilligen$refmat  <- refmat
uitdewilligen$sizemat <- sizemat
uitdewilligen$ploidy  <- ploidy
msub <- mupdog(refmat = refmat, sizemat = sizemat, ploidy = ploidy, verbose = FALSE, control = list(obj_tol = 10^-6))
mupout <- msub

save(uitdewilligen, file = "../data/uitdewilligen.RData")
save(mupout, file = "../data/uitdewilligen_fit.RData")


index <- 71
plot(msub, index)
plot(mout, which_snp[index])

plot(msub$bias, mout$bias[which_snp])
cor(msub$bias, mout$bias[which_snp])
