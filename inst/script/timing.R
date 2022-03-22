load_all()
y <- makeSimSwishData(m=5000,n=20,null=TRUE)
y <- scaleInfReps(y)

y$allele <- factor(rep(1:2,times=10))
y$pair <- factor(rep(1:10,times=2))
y$pair2 <- factor(rep(1:10,each=2))
y$continuous <- round(rep(rnorm(10),each=2),2)

# 1.5 s
system.time({
  y <- swish(y, x="condition", quiet=TRUE)
})
# 26 s
system.time({
  y <- swish(y, x="condition", pair="pair", quiet=TRUE)
})
# 2.0 s
system.time({
  y <- swish(y, x="allele", pair="pair2", cov="condition", interaction=TRUE, quiet=TRUE)
})
# 3.0 s
system.time({
  y <- swish(y, x="allele", pair="pair2", cov="continuous", cor="pearson", quiet=TRUE)
})
