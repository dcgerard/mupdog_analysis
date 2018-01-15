## Format updog simulations for presentations better.
suppressMessages(library(tidyverse))
dat <- as_data_frame(read.csv("../../../reproduce_genotyping/Output/sims_out/sims_out.csv", row.names = NULL))


longdat <- dat %>% transmute(updog = uham, Blischak = bham, allele_freq = allele_freq, od_param = od_param, bias_val = bias_val) %>%
  gather(key = "Method", value = "PropCorrect", updog:Blischak)


pl <- ggplot(data = filter(longdat, bias_val == 1), mapping = aes(y = PropCorrect, x = allele_freq, color = Method, group = Method)) +
  facet_grid(. ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Proportion Correct") +
  geom_point(size = 0.1) +
  geom_smooth(color = "black") +
  ggtitle("Bias = 1")
pdf(file = "../output/updog_sims/bias1.pdf",
    family = "Times", colormodel = "cmyk", height = 3, width = 9.5)
print(pl)
dev.off()

pl <- ggplot(data = filter(longdat, bias_val == 0.75), mapping = aes(y = PropCorrect, x = allele_freq, color = Method, group = Method)) +
  facet_grid(. ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Proportion Correct") +
  geom_point(size = 0.1) +
  geom_smooth(color = "black") +
  ggtitle("Bias = 0.75")
pdf(file = "../output/updog_sims/bias75.pdf",
    family = "Times", colormodel = "cmyk", height = 3, width = 9.5)
print(pl)
dev.off()

pl <- ggplot(data = filter(longdat, bias_val == 0.5), mapping = aes(y = PropCorrect, x = allele_freq, color = Method, group = Method)) +
  facet_grid(. ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Proportion Correct") +
  geom_point(size = 0.1) +
  geom_smooth(color = "black") +
  ggtitle("Bias = 0.5")
pdf(file = "../output/updog_sims/bias5.pdf",
    family = "Times", colormodel = "cmyk", height = 3, width = 9.5)
print(pl)
dev.off()


longdat <- dat %>% select(od_param, uod_param, bias_val, ubias_val, seq_error, useq_error, allele_freq)



pl_od <- ggplot(data = longdat, mapping = aes(x = as.factor(od_param), y = uod_param)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(hat(tau))) +
  geom_hline(yintercept = 0, lty = 2, alpha = 1 / 2, color = "red") +
  geom_hline(yintercept = 0.01, lty = 2, alpha = 1 / 2, color = "red") +
  geom_hline(yintercept = 0.1, lty = 2, alpha = 1 / 2, color = "red") +
  ggtitle("Overdispersion Estimates")

pdf(file = "../output/updog_sims/od_est.pdf",
    family = "Times", colormodel = "cmyk", height = 3, width = 3)
print(pl_od)
dev.off()


pl_bias <- ggplot(data = longdat, mapping = aes(x = as.factor(od_param), y = log2(ubias_val))) +
  facet_grid(.~bias_val) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(log[2](hat(h)))) +
  geom_hline(mapping = aes(yintercept = log2(bias_val)), lty = 2, color = "red") +
  ggtitle("Bias Estimates")

pdf(file = "../output/updog_sims/bias_est.pdf",
    family = "Times", colormodel = "cmyk", height = 3, width = 3)
print(pl_bias)
dev.off()


pl_seq <- ggplot(data = longdat, mapping = aes(x = as.factor(od_param), y = useq_error)) +
  geom_boxplot(outlier.size = 0.1) +
  geom_hline(yintercept = 0.005, col = "red", lty = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(hat(epsilon))) +
  ggtitle("Sequencing Error Rate Estimates")

pdf(file = "../output/updog_sims/seq_est.pdf",
    family = "Times", colormodel = "cmyk", height = 3, width = 3)
print(pl_seq)
dev.off()


## Updog features

shir <- read.csv("../output/updog_sims/shir_features.csv")
pl <- ggplot(data = filter(shir, Parameter == "Log2-Bias"), mapping = aes(x = Value)) +
  geom_histogram(bins = 20, fill = "white", color = "black") +
  xlab("Log2-bias") +
  theme_bw()
pdf(file = "../output/updog_sims/shir_bias.pdf",
    family = "Times", colormodel = "cmyk", height = 3, width = 3)
print(pl)
dev.off()
