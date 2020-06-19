library(tidyverse)
library(RColorBrewer)

method = "lasso" # change to "boost" or "kernel" for plotting those

colors = brewer.pal(8, "Dark2")
if (method == "boost") {
  breaks = c("S", "T", "X", "U", "CB", "R")
  colors = c(colors[2], colors[1], colors[3], colors[7], colors[8], colors[4])
  shapes = c(11,15,12,13,17,8)
  sizes = c(1.5,2,2,2,2,1.7)
} else if (method == "lasso") {
  breaks = c("S", "T", "X", "U", "RS", "R")
  colors = c(colors[1], colors[5], colors[3], colors[7], colors[8], colors[4])
  shapes = c(15,18,12,13,17,8)
  sizes = c(2,2.8,2,2,2,1.7)
} else if (method == "kernel") {
  breaks = c("S", "T", "X", "U", "R")
  colors = c(colors[1], colors[5], colors[3], colors[7], colors[4])
  shapes = c(15,18,12,13,8)
  sizes = c(2,2.8,2,2,1.7)
}

plotsize = function(x,y) options(repr.plot.width=x, repr.plot.height=y)

plotsize(8,8)
out = read_csv(paste0("output_", method, ".csv"))


label_wrap <- function(variable, value) {
  paste0("Setup ", value)
}
out %>%
  select(-n, -d, -sigma, -X1) %>%
  gather(learner, mse, -setup, -oracle) %>%
  ggplot(aes(x=log(oracle), y=log(mse), color=learner, shape=learner, size=learner,breaks=learner)) +
  scale_size_manual(breaks=breaks, values=sizes)+
  scale_shape_manual(breaks=breaks, values=shapes)+
  scale_colour_manual(breaks=breaks, values = colors) +
  geom_point() +
  geom_abline(slope=1) +
  facet_wrap(~setup, ncol=2, scales="free", labeller = label_wrap) +
  labs(x = "oracle log mean-squared error", y = "log mean-squared error")
ggsave(paste0("plots/fig_", method, ".pdf"))

