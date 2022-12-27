library(dplyr)
library(tidyverse)
library(RColorBrewer)

#select the best BLAST hit based on e-value and bit score
setwd("/PathToProcessedAmpliconReads/")
allfiles = list.files(pattern = "*_unique.txt")

final <- c();
for (f in 1:length(allfiles)) {
  round <- read.table(allfiles[f], header = FALSE, sep = "\t")
  names(round) <- c("ID", "amplicon", "identity", "alignment.length", "mismatches", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")
  round$replicate <- substr(allfiles[f], start = 1, stop = 9)
  
  best.round <- c();
  for(i in 1:length(round[duplicated(round$ID), 1])){
    A <- subset(round, ID == round[duplicated(round$ID), 1][i])
    B <- subset(A, evalue == min(evalue) & max(bit.score))
    C <- length(B$ID) > 1
    
    if (C == FALSE){
      best.round <- rbind(best.round, B) 
    }
  }
  
  D <- round[!duplicated(round$ID), ]
  
  final <- rbind(final, D, best.round)
}


#read in deletion allele files
delfiles = list.files(pattern = "*_DEL.txt")

deletion <- c();
for (f in 1:length(delfiles)) {
  x <- read.table(delfiles[f], header = FALSE, sep = "\t")
  names(x) <- c("ID", "amplicon", "identity", "alignment.length", "mismatches", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")
  x$replicate <- substr(allfiles[f], start = 1, stop = 9)
  
  deletion <- rbind(deletion, x)
}


#merge
amplicons <- rbind(final, deletion)

amplicons$replicate2 <- substr(amplicons$replicate, start = 1, stop = 6)
amplicons$replicate2[amplicons$replicate2 == "F1_R1_"] <- "CTL1"
amplicons$replicate2[amplicons$replicate2 == "F1_R2_"] <- "CTL2"


#data checking
qc <- as.data.frame(summarise(group_by(amplicons, amplicon), length(identity), mean(identity), sd(identity), mean(alignment.length), sd(alignment.length), min(alignment.length),
                              mean(mismatches), sd(mismatches)))
#total reads: 2,268,810
#wild-type: 1,586,062
#deletion: 682,748


#tabulate allele counts and calculate frequencies
freq.table <- c();
for (i in 1:length(levels(as.factor(amplicons$replicate2)))) {
  A <- amplicons[amplicons$replicate2 == levels(as.factor(amplicons$replicate2))[i], ];
  
  B <- as.data.frame(table(A$amplicon))
  
  freq.table <- rbind(freq.table, data.frame(ID  = A$replicate2[1],
                                             DEL = B[1, 2],
                                             WT  = B[2, 2]))
}

freq.table$TOTAL <- rowSums(freq.table[, 2:3])
freq.table$F.DEL <- freq.table$DEL / freq.table$TOTAL
freq.table$F.WT  <- freq.table$WT / freq.table$TOTAL
freq.table$GEN   <- c(rep(0, 2), rep(11, 10), rep(20, 10))
freq.table$REP   <- substr(freq.table$ID, start = 4, stop = 6)

x <- c(0, 500, 500, 1000, 0.5, 0.5, 0, 0)

freq.table <- rbind(freq.table, x)

--------------------------------------------------------------------------------------------------------------

#estimate selection coefficient for locus
#Model 1: Estimate selection coefficient using 50% deletion allele at generation 1
selection1 <- glm(cbind(WT, DEL) ~ GEN, family = binomial(link = "logit"), data = freq.table[3:23, ])
summary(selection1)

#Model 2: Estimate selection coefficient using generations 11 and 20 only
selection2 <- glm(cbind(WT, DEL) ~ GEN, family = binomial(link = "logit"), data = freq.table[3:22, ])
summary(selection2)

--------------------------------------------------------------------------------------------------------------

# Estimate allele frequency trajectories under genetic drift
# Input variables:
# N     number of individuals of each sex: total population = 2N
# x.m   starting population of males by diploid genotype
# x.f   starting population of females by diploid genotype
# sm    cost of the A1 allele in males
# sf    cost of the A2 allele in females
# gens  number of generations

simulation <- function(N, x.m, x.f, gens){
  
  # record: generation, number of individuals of each diploid genotype by sex
  
  # selection vector by diploid genotype
  w.m <- as.vector(c(1, 1, 1))
  w.f <- as.vector(c(1, 1, 1))
  
  # number of individuals of each genotype after selection
  x.prime.m <- rbinom(3, x.m, w.m)
  x.prime.f <- rbinom(3, x.f, w.f)
  
  # frequency of each genotype after selection
  freq.x.prime.m <- matrix(x.prime.m / sum(x.prime.m), nrow = 3, ncol = 1)
  freq.x.prime.f <- matrix(x.prime.f / sum(x.prime.f), nrow = 1, ncol = 3)
  
  # frequency of each of the 9 possible matings
  mating.freq <- (freq.x.prime.m %*% freq.x.prime.f)
  
  # Determining each genotype frequency in zygotes: a11 = A1A1, a12 = A1A2, a22 = A2A2
  prob.a11 <- (mating.freq[1, 1] + 0.5 * mating.freq[1, 2] + 0.5 * mating.freq[2, 1] + 0.25 * mating.freq[2, 2]) /
    sum(mating.freq)
  prob.a12 <- (0.5 * mating.freq[1, 2] + mating.freq[1, 3] + 0.5 * mating.freq[2, 1] + 0.5 * mating.freq[2, 2] +
                 0.5 * mating.freq[2, 3] + mating.freq[3, 1] + 0.5 * mating.freq[3, 2]) / sum(mating.freq)
  prob.a22 <- (0.25 * mating.freq[2, 2] + 0.5 * mating.freq[2, 3] + 0.5 * mating.freq[3, 2] + mating.freq[3, 3]) /
    sum(mating.freq)
  
  probabilities <- cbind(prob.a11, prob.a12, prob.a22)
  
  # The number of male and female zygotes for each diploid genotype
  next.gen.m <- rmultinom(1, size = N, prob = probabilities)
  next.gen.f <- rmultinom(1, size = N, prob = probabilities)
  
  # The frequency of the A1 allele in each sex at the end of a generation
  freq.a1.m  <- (next.gen.m[1, 1] + 0.5 * next.gen.m[2,1]) / N
  freq.a1.f  <- (next.gen.f[1, 1] + 0.5 * next.gen.f[2,1]) / N
  
  results	<-	matrix(nrow = gens + 2, ncol = 9)
  results[1,]	<-	c(0, x.m, x.f, freq.a1.m, freq.a1.f)
  
  gens	<-	seq(1, gens + 1, by = 1)
  
  for (g in gens) {
    x.m	<-	next.gen.m
    x.f <-  next.gen.f
    
    x.prime.m <- rbinom(3, x.m, w.m)
    x.prime.f <- rbinom(3, x.f, w.f)
    
    freq.x.prime.m <- matrix(x.prime.m / sum(x.prime.m), nrow = 3, ncol = 1)
    freq.x.prime.f <- matrix(x.prime.f / sum(x.prime.f), nrow = 1, ncol = 3)
    
    mating.freq <- (freq.x.prime.m %*% freq.x.prime.f)
    
    prob.a11 <- (mating.freq[1, 1] + 0.5 * mating.freq[1, 2] + 0.5 * mating.freq[2, 1] +
                   0.25 * mating.freq[2, 2]) / sum(mating.freq)
    prob.a12 <- (0.5 * mating.freq[1, 2] + mating.freq[1, 3] + 0.5 * mating.freq[2, 1] + 0.5 * mating.freq[2, 2] +
                   0.5 * mating.freq[2, 3] + mating.freq[3, 1] + 0.5 * mating.freq[3, 2]) / sum(mating.freq)
    prob.a22 <- (0.25 * mating.freq[2, 2] + 0.5 * mating.freq[2, 3] + 0.5 * mating.freq[3, 2] +
                   mating.freq[3, 3]) / sum(mating.freq)
    
    probabilities <- cbind(prob.a11, prob.a12, prob.a22)
    
    next.gen.m <- rmultinom(1, size = N, prob = probabilities)
    next.gen.f <- rmultinom(1, size = N, prob = probabilities)
    
    freq.a1.m  <- (next.gen.m[1, 1] + 0.5 * next.gen.m[2,1]) / N
    freq.a1.f  <- (next.gen.f[1, 1] + 0.5 * next.gen.f[2,1]) / N
    
    results[g+1,]	<-	c(g, x.m[1, 1], x.m[2, 1], x.m[3, 1], x.f[1, 1], x.f[2, 1], x.f[3, 1], freq.a1.m, freq.a1.f)
    
  }
  results	<-	as.data.frame(results)
  names(results)	<-	c("gen", "A1A1.m", "A1A2.m", "A2A2.m", "A1A1.f","A1A2.f", "A2A2.f",
                      "Freq A1 Males", "Freq A1 Females")
  return(results)
}
