library(tidyverse)
plotsize = function(x,y) options(repr.plot.width=x, repr.plot.height=y)

plotsize(8,8)
out = read_csv("output_boost.csv")
out %>%
  select(-n, -p, -sigma, -X1) %>%
  gather(learner, mse, -setup, -oracle) %>%
  ggplot(aes(x=log(oracle), y=log(mse), color=learner)) +
  geom_point() +
  geom_abline(slope=1) +
  facet_wrap(~setup, ncol=2, scales="free")
