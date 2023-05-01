# Load in raw cAMP data

library(ggplot2)
library(tidyverse)

cAMP <- read.csv("PATHtoSUM.csv")

cAMP$batch <- rownames(cAMP)

colnames(cAMP)[c(1,2,3)] <- c("WT", "pilU", "pilT")
colnames(cAMP)



cAMPlong <- gather(cAMP, key="Geno", value="Signal", -batch)

cAMPlong$Geno <- factor(cAMPlong$Geno)

cAMPlong$Geno <- relevel(cAMPlong$Geno, ref = "WT") 


# plot the data to see what's going on

# plot the data to see what's going on
ggplot(cAMPlong)+
  theme_dark()+
  scale_fill_viridis_d()+
  scale_color_viridis_d(option = "magma")+
  geom_violin(aes(y=Signal, x=Geno, group=Geno, fill=Geno))+
  geom_boxplot(width=0.1, aes(y=Signal, x=Geno, group=Geno, col=batch)) +
  geom_point(size = 3, aes(y=Signal, x=Geno, group=Geno, col=batch)) +
  theme(axis.title.x = element_text(size=15))+
  theme(axis.title.y = element_text(size=15))+
  theme(axis.text.y = element_text(size=12, face = "bold"))+
  theme(axis.text.x = element_text(size=12, face = "bold"))+
  ggtitle("cAMP by genotype and batch")


# Calculate differences from WT
cAMPmodel <- lm(Signal ~ Geno + batch, data = cAMPlong)

summary(cAMPmodel)

cAMPmodelCoeff <- data.frame(coefficients(summary(cAMPmodel)))

# Load in Twitch Diameter

TwitchD <- read.delim("PATHtoSUM.txt", sep = "\t")
colnames(TwitchD)

TwitchD$batch <- rownames(TwitchD)

TwitchDlong <- gather(TwitchD, key="Geno", value="D", -batch)

TwitchDlong$Geno <- factor(TwitchDlong$Geno)

TwitchDlong$Geno <- relevel(TwitchDlong$Geno, ref = "WT")


# plot the data to see what's going on
ggplot(TwitchDlong)+
  theme_dark()+
  scale_fill_viridis_d()+
  scale_color_viridis_d(option = "magma")+
  geom_violin(aes(y=D, x=Geno, group=Geno, fill=Geno))+
  geom_boxplot(width=0.1, aes(y=D, x=Geno, group=Geno, col=batch)) +
  geom_point(size = 3, aes(y=D, x=Geno, group=Geno, col=batch)) +
  theme(axis.title.x = element_text(size=15))+
  theme(axis.title.y = element_text(size=15))+
  theme(axis.text.y = element_text(size=12, face = "bold"))+
  theme(axis.text.x = element_text(size=12, face = "bold"))+
  ggtitle("Twitch Diameter by genotype and batch")

# Make a model

TwitchDmodel <- lm(D ~ Geno + batch, data = TwitchDlong)

summary(TwitchDmodel)

TwitchDmodelCoeff <- data.frame(coefficients(summary(TwitchDmodel)))

# Interaction data

Ixn <- read.csv("PATHtoSUM.csv")

Ixn$batch <- rownames(Ixn)

colnames(Ixn)[c(1)] <- c("WT")

Ixnlong <- gather(Ixn, key="Geno", value="Ixn", -batch)

Ixnlong$Geno <- factor(Ixnlong$Geno)

Ixnlong$Geno <- relevel(Ixnlong$Geno, ref = "WT")

ggplot(Ixnlong)+
  theme_dark()+
  scale_fill_viridis_d()+
  scale_color_viridis_d(option = "magma")+
  geom_violin(aes(y=Ixn, x=Geno, group=Geno, fill=Geno))+
  geom_boxplot(width=0.1, aes(y=Ixn, x=Geno, group=Geno)) +
  geom_point(size = 3, aes(y=Ixn, x=Geno, group=Geno, col=batch)) +
  theme(axis.title.x = element_text(size=15))+
  theme(axis.title.y = element_text(size=15))+
  theme(axis.text.y = element_text(size=12, face = "bold"))+
  theme(axis.text.x = element_text(size=12, face = "bold"))+
  ggtitle("PilT-PilJ Ixn by genotype and batch")

# Make a model


IxnModel <- lm(Ixn ~ Geno + batch, data = Ixnlong)

summary(IxnModel)

IxnmodelCoeff <- data.frame(coefficients(summary(IxnModel)))

# Combine all three models

cAMP.Twitch <- merge(cAMPmodelCoeff[1:8,], TwitchDmodelCoeff[1:8,], by = 0)

colnames(cAMP.Twitch) <- sub("*.x$", ".cAMP",  colnames(cAMP.Twitch))
colnames(cAMP.Twitch) <- sub("*.y$", ".twitch",  colnames(cAMP.Twitch))

cAMP.Twitch.Ixn <- merge(cAMP.Twitch, IxnmodelCoeff[1:8,], by.x = "Row.names", by.y = 0 )

colnames(cAMP.Twitch.Ixn )[10:13 ] <- paste(colnames(cAMP.Twitch.Ixn )[10:13 ], "Ixn", sep = '.')

colnames(cAMP.Twitch.Ixn )[1] <- "Genotype"

cAMP.Twitch.Ixn$Genotype[1] <- "WT"
cAMP.Twitch.Ixn$Genotype[2] <- "D31K"
cAMP.Twitch.Ixn$Genotype[3] <- "E204A"
cAMP.Twitch.Ixn$Genotype[4] <- "H222A"
cAMP.Twitch.Ixn$Genotype[5] <- "H229A"
cAMP.Twitch.Ixn$Genotype[6] <- "H44L"
cAMP.Twitch.Ixn$Genotype[7] <- "K136A"
cAMP.Twitch.Ixn$Genotype[8] <- "K58A"


# All estimates are relative to WT. Add WT value to the others

cAMP.Twitch.Ixn$Estimate.cAMP[2:8] <- cAMP.Twitch.Ixn$Estimate.cAMP[2:8] + cAMP.Twitch.Ixn$Estimate.cAMP [1]
cAMP.Twitch.Ixn$Estimate.twitch [2:8] <- cAMP.Twitch.Ixn$Estimate.twitch[2:8] + cAMP.Twitch.Ixn$Estimate.twitch [1]
cAMP.Twitch.Ixn$Estimate.Ixn [2:8] <- cAMP.Twitch.Ixn$Estimate.Ixn [2:8] + cAMP.Twitch.Ixn$Estimate.Ixn [1]

plot(Estimate.cAMP ~ Estimate.twitch, data = cAMP.Twitch.Ixn)
abline(lm(Estimate.cAMP ~ Estimate.twitch, data = cAMP.Twitch.Ixn))

cAMP.Twitch.Ixn$Genotype <- factor(cAMP.Twitch.Ixn$Genotype)

#cAMP and twitch
ggplot(cAMP.Twitch.Ixn, aes(y = Estimate.cAMP, x = Estimate.twitch, label=cAMP.Twitch.Ixn$Genotype) ) + geom_point() + 
  geom_point(size = 4) + 
  geom_smooth(method = 'lm') + 
  geom_text(size=5, nudge_y = 0, vjust=-1) +
  theme(axis.title.x = element_text(size=15))+
  theme(axis.title.y = element_text(size=15))+
  theme(axis.text.y = element_text(size=12, face = "bold"))+
  theme(axis.text.x = element_text(size=12, face = "bold"))+
  ggtitle("cAMP and Twitch Diameter: Rsq:.9579, AdjRsq:.9509, pval:2.369e-5")
ggsave("PATHtoSUM.pdf")

summary(lm(Estimate.cAMP ~ Estimate.twitch, data = cAMP.Twitch.Ixn))



#cAMP and T-J interaction
plot(Estimate.cAMP ~  Estimate.Ixn, data = cAMP.Twitch.Ixn)
abline(lm(Estimate.cAMP ~  Estimate.Ixn, data = cAMP.Twitch.Ixn))

ggplot(cAMP.Twitch.Ixn, aes(y = Estimate.cAMP, x = Estimate.Ixn, label=cAMP.Twitch.Ixn$Genotype) ) + geom_point() + 
  geom_point(size = 4) + 
  geom_smooth(method = 'lm') + 
  geom_text(size=5, nudge_y = 0, vjust=-1) +
  theme(axis.title.x = element_text(size=15))+
  theme(axis.title.y = element_text(size=15))+
  theme(axis.text.y = element_text(size=12, face = "bold"))+
  theme(axis.text.x = element_text(size=12, face = "bold"))+
  ggtitle("cAMP and PilT-PilJ Interaction: Rsq: .228, AdjRsq: .09938, pval: .2314")

ggsave("PATHtoSUM.pdf")

summary(lm(Estimate.cAMP ~ Estimate.Ixn, data = cAMP.Twitch.Ixn))

