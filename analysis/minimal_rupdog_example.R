library(updog)
nind              <- 1000

uobj              <- list()
class(uobj)       <- "updog"
uobj$input        <- list()
uobj$input$model  <- "f1"
uobj$input$ploidy <- 6
uobj$input$osize  <- rpois(n = nind, lambda = 100)
uobj$input$p1size <- 100
uobj$input$p2size <- 100
uobj$p1geno       <- 3
uobj$p2geno       <- 3
uobj$bias_val     <- 1
uobj$od_param     <- 0.001
uobj$seq_error    <- 0.005
uobj$allele_freq  <- -1 ## hacky thing I need but should get rid of
uobj$out_prop     <- 0
uobj$out_mean     <- 1/2
uobj$out_disp     <- 1/3

rout <- rupdog(uobj)
plot(rout)
