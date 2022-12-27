library(dplyr)
library(tidyverse)

#Analyze Late Life Fecundity
fecundity <- read.table("/Pathto/FileS1_FecundityData.txt", header = TRUE, sep = "\t)

fecundity$strain <- c(rep("JK", 30), rep("VX", 30))
fecundity$total <- rowSums(fecundity[, 2:9], na.rm = F)

t.test(fecundity$total ~ fecundity$strain, alternative = "two.sided", na.action = na.omit)

power <- pwr.t2n.test(n1 = length(subset(fecundity, strain=="JK" & is.na(total)==FALSE)$total), n2 = length(subset(fecundity, strain=="VX" & is.na(total)==FALSE)$total),
                      d = NULL, sig.level = 0.05, power = 0.8)

----------------------------------------------------------------------------------------------------------------------------------------------------------

#Analyze Immunostaining Data
#Read in data files
setwd("/PathToFileS4_IHSdata.xlsx/")

fulldata <- 
  list.files(pattern = "*.txt") %>% 
  map_df(~read.table(., header = TRUE, sep = "\t"))

matedFemale <- subset(fulldata, Worm == "JK574(female)_VX300(male)_D2", c("Worm", "Head", "UpperGonad", "MidGonad", "LowerGonad", "VulvaCloaca", "Tail", "None", "Eggs"))
onlyHerm    <- subset(fulldata, Worm == "VX300_D1_herm_nomales", c("Worm", "Head", "UpperGonad", "MidGonad", "LowerGonad", "VulvaCloaca", "Tail", "None"))
data        <- subset(fulldata, Worm != "JK574(female)_VX300(male)_D2" & Worm != "VX300_D1_herm_nomales", 
                      c("Worm", "Head", "UpperGonad", "MidGonad", "LowerGonad", "VulvaCloaca", "Tail", "None", "Censor"))

data$strain <- substr(data$Worm, start = 1, stop = 5)
data$age    <- substr(data$Worm, start = 7, stop = 8)
data$sex    <- substr(data$Worm, start = 10, stop = 13)

data$sex[data$sex == "Herm"]   <- "herm"
data$sex[data$sex == "Male"]   <- "male"
data$sex[data$sex == ""]       <- "unknown"
data$Censor[is.na(data$Censor) == TRUE] <- 0

data.filtered <- subset(data, Censor != 1)

allcounts <- data.filtered %>%
              group_by(strain, sex, age) %>%
              summarise(length(Worm))
allcounts <- as.data.frame(allcounts)
names(allcounts) <- c("strain", "sex", "age", "total.individ")


#summarize female data
herm <- subset(data.filtered, sex == "herm") %>%
          group_by(strain, age) %>%
          summarise(sum(Head), sum(UpperGonad), sum(MidGonad), sum(LowerGonad), sum(VulvaCloaca), sum(Tail), sum(None))
herm <- as.data.frame(herm)

names(herm) <- c("strain", "age", "head", "upperG", "midG", "lowerG", "vulva", "tail", "none")

herm$total <- rowSums(herm[, 3:9])

l4.test.herm <- chisq.test(herm[herm$age == "L4", c(7, 9)])
d1.test.herm <- chisq.test(herm[herm$age == "D1", c(7, 9)])
d2.test.herm <- chisq.test(herm[herm$age == "D2", c(7, 9)])
d8.test.herm <- chisq.test(herm[herm$age == "D8", c(3:9)])


#summarize male data
male <- subset(data.filtered, sex == "male") %>%
  group_by(strain, age) %>%
  summarise(sum(Head), sum(UpperGonad), sum(MidGonad), sum(LowerGonad), sum(VulvaCloaca), sum(Tail), sum(None))
male <- as.data.frame(male)

names(male) <- c("strain", "age", "head", "upperG", "midG", "lowerG", "cloaca", "tail", "none")

male$total <- rowSums(male[, 3:9])

l4.test.male <- chisq.test(male[male$age == "L4", c(4, 5, 8, 9)])
d1.test.male <- chisq.test(male[male$age == "D1", c(8, 9)])
d2.test.male <- chisq.test(male[male$age == "D2", c(8, 9)])
d8.test.male <- chisq.test(male[male$age == "D8", c(8:9)])


#JK574 females x VX300 males
cross <- as.data.frame(summarise(group_by(matedFemale, Eggs), sum(Head), sum(UpperGonad), sum(MidGonad), sum(LowerGonad), sum(VulvaCloaca), sum(Tail), sum(None)))
names(cross) <- c("strain", "head", "upperG", "midG", "lowerG", "vulva", "tail", "none")

cross$total <- rowSums(cross[, 2:8])

cross <- rbind(cross, herm[herm$age == "D2", c(1, 3:10)])

cross.test1 <- chisq.test(cross[1:3, c(6, 8)])
cross.test2 <- chisq.test(cross[c(1, 3), c(6, 8)])
cross.test3 <- chisq.test(cross[2:3, c(6, 8)])


#VX300 Herms only
nomales <- as.data.frame(summarise(group_by(onlyHerm, Worm), sum(Head), sum(UpperGonad), sum(MidGonad), sum(LowerGonad), sum(VulvaCloaca), sum(Tail), sum(None)))
names(nomales) <- c("strain", "head", "upperG", "midG", "lowerG", "vulva", "tail", "none")

nomales$total <- rowSums(nomales[, 2:8])

nomales <- rbind(nomales, herm[herm$age == "D1", c(1, 3:10)])

nomales.test1 <- chisq.test(nomales[1:2, c(6, 8)])
nomales.test2 <- chisq.test(nomales[c(1, 3), c(6, 8)])
