library(tidyverse)

args=(commandArgs(TRUE))
method = as.character(args[1])

if (method == "boost") {
  breaks = c("S", "T", "X", "U", "CB", "R")
} else if (method == "lasso") {
  breaks = c("S", "T", "X", "U", "R", "RS")
}

plotsize = function(x,y) options(repr.plot.width=x, repr.plot.height=y)

plotsize(8,8)
out = read_csv(paste0("output_", method, ".csv"))
out %>%
  select(-n, -p, -sigma, -X1) %>%
  gather(learner, mse, -setup, -oracle) %>%
  ggplot(aes(x=log(oracle), y=log(mse), color=learner, shape=learner, breaks=learner)) +
  scale_shape_manual(breaks=breaks,
                     values=c(17,9,14,15,16,8))+
  scale_colour_manual(breaks=breaks,
                      values = c("purple","blue", "green", "red", "cyan", "magenta"))+
  geom_point() +
  geom_abline(slope=1) +
  facet_wrap(~setup, ncol=2, scales="free") +
  labs(x = "oracle log mean-squared error", y = "log mean-squared error")
ggsave(paste0("plots/fig_", method, ".pdf"))

