library(ggplot2)
library(gganimate)
library(SummarizedExperiment)
library(fishpond)
set.seed(5)
y <- makeSimSwishData(m=10, meanVariance=TRUE)
y$condition <- 1:10
assay(y, "mean")[1,] <- 100 + y$condition * 20 + rnorm(10,0,20)
assay(y, "variance")[1,] <- 200 + c(1,0,1,0,0, 1,0,0,0,1) * 10000
# plotInfReps(y, 1, x="condition")
y <- makeInfReps(y, numReps=30)
dat <- getTrace(y, idx=1, samp_idx=1:10)
library(dplyr)
coldata <- colData(y) %>%
  as.data.frame() %>%
  mutate(sample=condition)
dat <- dat %>% inner_join(coldata)
anim <- dat %>%
  ggplot(aes(condition, count)) +
  geom_point(col="dodgerblue", size=3) +
  geom_line(stat="smooth", method="lm", formula=y~x, se=FALSE,
            col="black", lwd=2) +
  ylab("scaled counts") +
  coord_cartesian(ylim=c(0,500)) + 
  theme_bw() +
  theme(text = element_text(size = 24)) +
  transition_time(infRep) +
  shadow_mark(alpha=0.25, size=2.5)
animate(anim, duration=6, fps=5)

library(tidyr)
library(purrr)
library(broom)
dat %>%
  select(-sample) %>%
  nest(data=-infRep) %>%
  mutate(fit = map(data, ~ tidy(lm(count ~ condition, data=.x)))) %>%
  unnest(fit) %>%
  group_by(term) %>%
  summarize(mu=mean(estimate))

metadata(y)$infRepsScaled <- TRUE
library(rafalib)
bigpar()
plotInfReps(y, idx=1, x="condition", main="")
abline(120, 17.7, col=rgb(0,0,0,.3), lwd=5)
