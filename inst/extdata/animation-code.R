library(ggplot2)
library(gganimate)
library(SummarizedExperiment)
library(fishpond)
set.seed(5)
y <- makeSimSwishData(m=10, meanVariance=TRUE)
y$condition <- 0:9
assay(y, "mean")[1,] <- 100 + y$condition * 20 + rnorm(10,0,20)
assay(y, "variance")[1,] <- 200 + c(0,1,0,0,0, 1,0,0,0,1) * 1000
plotInfReps(y, 1, x="condition")
y <- makeInfReps(y, numReps=20)
dat <- getTrace(y, idx=1, samp_idx=1:10)
library(dplyr)
coldata <- colData(y) %>%
  as.data.frame() %>%
  mutate(sample=1:10)
dat <- dat %>% inner_join(coldata)
dat %>% ggplot(aes(condition, count)) +
  geom_point(col="dodgerblue") +
#  geom_smooth(method="lm", se=FALSE) +
  transition_time(infRep)
