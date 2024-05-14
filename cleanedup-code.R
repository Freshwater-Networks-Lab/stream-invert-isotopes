#cleaned up code
library(rstatix)
library(ggtext)
library(ggplot2)
library(rjags)
library(SIBER)
library(dplyr)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(reshape2)
library(pals)
library(forcats)

#community analysis

#create stacked bar plots for each stream for each family and FFG
#2 plots for raw abundace, 2 for percentage representation

#FFG
guilds <- read.csv(file.choose())
guilds
guilds <- guilds[,-1]
guilds

positions <- c("Ll6", "Ll7", "Ll3", "Ll8")
ggplot(guilds, aes(x = variable, y = value, fill = guild))+
  scale_x_discrete(limits = positions) +
  labs(x = "Stream",
       y = "Abundance")+
  geom_col()

ggplot(guilds, aes(x = variable, y = value, fill = guild))+
  scale_x_discrete(limits = positions) +
  labs(x = "Stream",
       y = "Percentage")+
  geom_col(position ="fill")+
  scale_y_continuous(labels = scales::percent)

#Family
family <- read.csv(file.choose())
family2 <- melt(family, id = "Family")
family2
Family <- factor(family2$Family)
family2$Family <- factor(family2$Family, levels = c("Baetidae", "Heptageniidae", "Rhyacophalidae", "Polycentropodidae",
                                                    "Hydropsychidae", "Limnephilidae", "Philoptamidae", "Psychomyiidae", 
                                                    "Goeridae", "Perlodidae", "Chloroperlidae", "Leuctridae",
                                                    "Nemouridae", "Taeniopterygidae", "Simuliidae", "Pedicidae",
                                                    "Scirtridae", "Elmidae", "Hydraenidae", "Chironomidae", "Valvatidae") )
#raw abundance

ggplot(family2, aes(x = variable, y = value, fill = Family)) +
  scale_fill_manual(values = unname(alphabet()))+
  labs(x="Stream",
       y="Abundance") +
  geom_col()

#percentage
ggplot(family2, aes(x = variable, y = value, fill = Family)) +
  scale_fill_manual(values = unname(alphabet()))+
  labs(x="Stream",
       y="Abundance") +
  geom_col(position ="fill")+
  scale_y_continuous(labels = scales::percent)


#wilcoxon test for community data
familystats <- family[,2:5]
familystats

wilcox.test(familystats$Ll3, familystats$Ll6, exact = FALSE)
wilcox.test(familystats$Ll3, familystats$Ll7, exact = FALSE)
wilcox.test(familystats$Ll3, familystats$Ll8, exact = FALSE)
wilcox.test(familystats$Ll6, familystats$Ll7, exact = FALSE)
wilcox.test(familystats$Ll6, familystats$Ll8, exact = FALSE)
wilcox.test(familystats$Ll7, familystats$Ll8, exact = FALSE)









#####
#stable isotope analysis

#GLM
#read in file all-data
overall <- read.csv(file.choose())
#look at data
overall
#get rid of NAs and then have a look
overall <- na.omit(overall)
overall

#GLM for stream site
model1 <- glm(iso1 ~ community, data = overall, family = gaussian (link = "log"))
car::Anova(model1, type=3)
summary.lm(model1)
model2 <- glm(-iso2 ~ community, data = overall, family = gaussian (link = "log"))
car::Anova(model2, type=3)
summary.lm(model2)

#GLM for FFG
model3 <- glm(iso1 ~ guild, data = overall, family = gaussian (link = "log"))
summary.lm(model3)
car::Anova(model3, type=3)

model4 <- glm(-iso2 ~ guild, data = overall, family = gaussian (link = "log"))
summary.lm(model4)
car::Anova(model4, type=3)



#####
#8 box plots - 2 per stream for carbon and nitrogen 

library(ggplot2)
#Ll3 nitrogen stream
new_order1 <- with(Ll3, reorder(group, iso1, median , na.rm=T))
new_order2 <- with(Ll3, reorder(group, iso2, median , na.rm=T))
ggplot(Ll3, aes(x=new_order1, y=iso1), )+
  labs(x="Taxon",
       y="δ<sup>15</sup>N",
       title = "Ll3") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y=element_markdown())+
  geom_boxplot()

ggplot(Ll3, aes(x=new_order2, y=iso2), )+
  labs(x="Taxon",
       y="δ<sup>13</sup>C",
       title = "Ll3") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_markdown())+
  geom_boxplot()

#Ll8
new_order8n <- with(Ll8clean2, reorder(group, iso1, median , na.rm=T))
new_order8c <- with(Ll8clean2, reorder(group, iso2, median , na.rm=T))

ggplot(Ll8clean2, aes(x=new_order8n, y=iso1), )+
  labs(x="Taxon",
       y="δ<sup>15</sup>N",
       title = "Ll8") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y=element_markdown())+
  geom_boxplot()

ggplot(Ll8clean2, aes(x=new_order8c, y=iso2), )+
  labs(x="Taxon",
       y="δ<sup>13</sup>C",
       title = "Ll8") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_markdown())+
  geom_boxplot()


#Ll6
new_order6n <- with(Ll62, reorder(group, iso1, median , na.rm=T))
new_order6c <- with(Ll62, reorder(group, iso2, median , na.rm=T))

ggplot(Ll62, aes(x=new_order6n, y=iso1), )+
  labs(x="Taxon",
       y="δ<sup>15</sup>N",
       title = "Ll6") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y=element_markdown())+
  geom_boxplot()

ggplot(Ll62, aes(x=new_order6c, y=iso2), )+
  labs(x="Taxon",
       y="δ<sup>13</sup>C",
       title = "Ll6") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_markdown())+
  geom_boxplot()



#Ll7
new_order7n <- with(Ll72, reorder(group, iso1, median , na.rm=T))
new_order7c <- with(Ll72, reorder(group, iso2, median , na.rm=T))

ggplot(Ll72, aes(x=new_order7n, y=iso1), )+
  labs(x="Taxon",
       y="δ<sup>15</sup>N",
       title = "Ll7") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.y=element_markdown())+
  geom_boxplot()

ggplot(Ll72, aes(x=new_order7c, y=iso2), )+
  labs(x="Taxon",
       y="δ<sup>13</sup>C",
       title = "Ll7") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y=element_markdown())+
  geom_boxplot()


#####
#intra analysis
#read in intra-data.csv
intradata <- read.csv(file.choose())
intradata
Ll6intra <- subset(intradata, community == "Ll6")
Ll7intra <-subset(intradata, community == "Ll7")
Ll3intra <-subset(intradata, community == "Ll3")
Ll8intra <-subset(intradata, community == "Ll8")
Ll8intra

#overall plot
object <- createSiberObject(intradata)

#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                main = "Ll6",
                x.limits = c(-34, -20),
                y.limits = c(0, 9))
object$group.names
a.name.order <- c("Baetidae", "Heptageniidae", "Leuctridae", "Protonemura")
legend("topright", a.name.order, 
       pch = c(1,1), col = c(1,2), lty = 1)
#Ll6 plot
object <- createSiberObject(Ll6intra)

#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                main = "Ll6",
                x.limits = c(-32, -21),
                y.limits = c(0, 6))
object$group.names
a.name.order <- c("Baetidae", "Heptageniidae")
legend("topright", a.name.order, 
       pch = c(1,1), col = c(1,2), lty = 1)


#Ll7 plot
object <- createSiberObject(Ll7intra)

#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                main = "Ll7",
                x.limits = c(-32, -21),
                y.limits = c(0, 6))
object$group.names
a.name.order <- c("Baetidae", "Heptageniidae")
legend("topright", a.name.order, 
       pch = c(1,1), col = c(1,2), lty = 1)

#Ll3 plot
object <- createSiberObject(Ll3intra)
palette(viridis(3))
#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                main = "Ll3",
                x.limits = c(-32, -24),
                y.limits = c(2, 7))
object$group.names
a.name.order <- c("Leuctridae")
legend("topright", a.name.order, 
       pch = c(1,1), col = c(1,2), lty = 1)
#Ll8 plot
object <- createSiberObject(Ll8intra)

#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                main = "Ll8",
                x.limits = c(-34, -24),
                y.limits = c(3, 9))

object$group.names
a.name.order <- c("Leuctridae", "Protonemura")
legend("topright", a.name.order, 
       pch = c(1,1), col = c(1,2), lty = 1)

library(SIBER)


#area intra
#define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.example, parms, priors)


ellipses.posterior <- siberMVN(object, parms, priors)

SEA.B <- siberEllipses(ellipses.posterior)
SEA.B

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('permille' ^2) ),
                 ylims = c(0,5),
                 bty = "L",
                 las = 1,
)



#wilcox test for ellipses sizes
#
SEA.B1 <- SEA.B[,1:2]
leuctridaell3 <- SEA.B[,1]
heptall6 <- SEA.B[,2]
baetidaell6 <- SEA.B[,3]
heptall7 <- SEA.B[,4]
baetidaell7 <- SEA.B[,5]
leuctridaell8 <- SEA.B[,6]
protoll8 <- SEA.B[,7]

wilcox.test(leuctridaell3, heptall6)
wilcox.test(leuctridaell3, baetidaell6)
wilcox.test(leuctridaell3, heptall7)
wilcox.test(leuctridaell3, baetidaell7)
wilcox.test(leuctridaell3, leuctridaell8)
wilcox.test(leuctridaell3, protoll8)

wilcox.test(heptall6, baetidaell6)
wilcox.test(heptall6, heptall7)
wilcox.test(heptall6, baetidaell7)
wilcox.test(heptall6, leuctridaell8)
wilcox.test(heptall6, protoll8)

wilcox.test(baetidaell6, heptall7)
wilcox.test(baetidaell6, baetidaell7)
wilcox.test(baetidaell6, leuctridaell8)
wilcox.test(baetidaell6, protoll8)

wilcox.test(heptall7, baetidaell7)
wilcox.test(heptall7, leuctridaell8)
wilcox.test(heptall7, protoll8)

wilcox.test(baetidaell7, leuctridaell8)
wilcox.test(baetidaell7, protoll8)

wilcox.test(leuctridaell8, protoll8)

summary(SEA.B)

#euclidean distance between every datapoint per taxon
#
Ll6intra <- subset(intradata, community == "Ll6")
Ll7intra <-subset(intradata, community == "Ll7")
Ll3intra <-subset(intradata, community == "Ll3")
Ll8intra <-subset(intradata, community == "Ll8")

Ll6intrah <- subset(Ll6intra, group == "heptageniidae")
Ll6comph <- Ll6intrah[,1:2]
Ll6intrab <- subset(Ll6intra, group == "baetidae")
Ll6compb <- Ll6intrab[,1:2]

Ll7intrah <- subset(Ll7intra, group == "heptageniidae")
Ll7comph <- Ll7intrah[,1:2]
Ll7intrab <- subset(Ll7intra, group == "baetidae")
Ll7compb <- Ll7intrab[,1:2]

Ll8intral <- subset(Ll8intra, group == "leuctridae")
Ll8compl <- Ll8intra[,1:2]
Ll8intrap <- subset(Ll8intra, group == "protonemura")
Ll8compp <- Ll8intrap[,1:2]
library(rdist)

Ll3distances <- rdist(Ll3comp)

Ll6distancesh <-rdist(Ll6comph)                     
Ll6distancesb <- rdist(Ll6compb) 

Ll7distancesh <-rdist(Ll7comph)                     
Ll7distancesb <- rdist(Ll7compb) 

Ll8distancesl <-rdist(Ll8compl)                     
Ll8distancesp <- rdist(Ll8compp) 

#now doing t.tests between ones in the same stream and between same taxa

#Ll3 leuctridae to Ll8 leuctridae

shapiro.test(intradata$iso2)

wilcox.test(Ll3distances, Ll8distancesl)

#data:  Ll3distances and Ll8distancesl
#W = 1741, p-value = 0.0209
#alternative hypothesis: true location shift is not equal to 0

#compare Ll8
wilcox.test(Ll8distancesl, Ll8distancesp)
#W = 3196, p-value = 0.7191
#alternative hypothesis: true location shift is not equal to 0

#compare Ll6
wilcox.test(Ll6distancesb, Ll6distancesh)
#W = 1025, p-value = 0.9232
#alternative hypothesis: true location shift is not equal to 0

#compare Ll7

wilcox.test(Ll7distancesh, Ll7distancesb)
#W = 820, p-value = 0.1217
#alternative hypothesis: true location shift is not equal to 0

#compare heptageniidae

wilcox.test(Ll6distancesh, Ll7distancesh)
#W = 963, p-value = 0.6939
#alternative hypothesis: true location shift is not equal to 0

#compare baetidae

wilcox.test(Ll6distancesb, Ll7distancesb)
#W = 744, p-value = 0.03009
#alternative hypothesis: true location shift is not equal to 0



#means of distances
summary(Ll3distances)
#1.4813
summary(Ll6distancesb)
#1.4737
summary(Ll6distancesh)
#1.5198
summary(Ll7distancesb)
#2.143
summary(Ll7distancesh)
# 1.7003
summary(Ll8distancesl)
#1.967
summary(Ll8distancesp)
#1.8334





#inter-specific analyses

#separate centroid and siber biplot by stream
#read all-siber
all3 <- read.csv(file.choose())
all3
Ll3 <- subset(all3, community == "Ll3")
Ll3
library(SIBER)
palette(viridis(5))
#siber biplot Ll3
object <-createSiberObject(Ll3)
#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                x.limits = c(-33, -21.5),
                y.limits = c(0, 9))

object$group.names
a.name.order <- c("Chloroperlidae", "Leuctridae", "Limnephilidae", "Perlodidae", "Polycentropodidae")
legend("topright", a.name.order, 
       pch = c(1,1,1,1,1), col = c(1,2,3,4,5), lty = 1)


#centroid plot
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(rstatix)
library(GGally)
library(plyr)

data.sum <- ddply(Ll3, c("group"), summarise,
                  iso1mn=mean(iso1),
                  iso1sd=sd(iso1),
                  iso2mn=mean(iso2),
                  iso2sd=sd(iso2))

data.sum


xlims <- aes(xmax = iso2mn + iso2sd, xmin = iso2mn - iso2sd)
ylims <- aes(ymax = iso1mn + iso1sd, ymin = iso1mn - iso1sd)

sum.biplot <- ggplot(data.sum, aes(x=iso2mn, y=iso1mn, colour = group)) +
  geom_point(size=3)+
  geom_errorbar(ylims, width=0.2)+
  geom_errorbarh(xlims, height=0.2) +
  ylab(expression(delta^{15}~N)) +
  xlab(expression(delta^{13}~C)) +
  scale_color_manual(values = c("chloroperlidae" = "black",
                                "leuctridae" = "red",
                                "limnephilidae" = "green",
                                "perlodidae" = "darkblue",
                                "polycentropodidae" = "cyan"
  )) +
  labs(title = "Ll3") +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title =element_text(size = 12),
        legend.text = element_text(size = 10.5),
        axis.line = element_line(colour = "black", linewidth = 0.6))
sum.biplot


#####
#Ll8
Ll8 <-subset(all3, community == "Ll8")
Ll8
Ll8clean2 <- Ll8[-c(23:24), ]
Ll8clean2
library(SIBER)
library(viridis)
palette(viridis(3))
#siber biplot Ll8
object <-createSiberObject(Ll8clean2)
#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                x.limits = c(-35, -24),
                y.limits = c(4, 9))

object$group.names
a.name.order <- c("Lecutridae", "Protonemura", "Tipulidae")
legend("topright", a.name.order, 
       pch = c(1,1,1), col = c(1,2,3), lty = 1)


#centroid plot
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(rstatix)
library(GGally)
library(plyr)

data.sum <- ddply(Ll8clean2, c("group"), summarise,
                  iso1mn=mean(iso1),
                  iso1sd=sd(iso1),
                  iso2mn=mean(iso2),
                  iso2sd=sd(iso2))

data.sum


xlims <- aes(xmax = iso2mn + iso2sd, xmin = iso2mn - iso2sd)
ylims <- aes(ymax = iso1mn + iso1sd, ymin = iso1mn - iso1sd)

sum.biplot <- ggplot(data.sum, aes(x=iso2mn, y=iso1mn, colour = group)) +
  geom_point(size=3)+
  geom_errorbar(ylims, width=0.2)+
  geom_errorbarh(xlims, height=0.2) +
  ylab(expression(delta^{15}~N)) +
  xlab(expression(delta^{13}~C)) +
  scale_color_manual(values = c("leuctridae" = "purple",
                                "protonemura" = "blue",
                                "tipulid" = "yellow"
  )) +
  labs(title = "Ll8") +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title =element_text(size = 12),
        legend.text = element_text(size = 10.5),
        axis.line = element_line(colour = "black", linewidth = 0.6))
sum.biplot

#####
#Ll6
Ll6 <-subset(all3, community == "Ll6")
Ll62 <- subset(Ll6, group == "heptageniidae" | group == "baetidae")
Ll62
#siber biplot Ll6
object <-createSiberObject(Ll62)
#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                x.limits = c(-32, -21),
                y.limits = c(1, 6))

object$group.names
a.name.order <- c("Baetidae", "Heptageniidae")
legend("topright", a.name.order, 
       pch = c(1,1), col = c(1,2), lty = 1)

#Ll6 centroid plot
data.sum <- ddply(Ll62, c("group"), summarise,
                  iso1mn=mean(iso1),
                  iso1sd=sd(iso1),
                  iso2mn=mean(iso2),
                  iso2sd=sd(iso2))

data.sum


xlims <- aes(xmax = iso2mn + iso2sd, xmin = iso2mn - iso2sd)
ylims <- aes(ymax = iso1mn + iso1sd, ymin = iso1mn - iso1sd)

sum.biplot <- ggplot(data.sum, aes(x=iso2mn, y=iso1mn, colour = group)) +
  geom_point(size=3)+
  geom_errorbar(ylims, width=0.2)+
  geom_errorbarh(xlims, height=0.2) +
  ylab(expression(delta^{15}~N)) +
  xlab(expression(delta^{13}~C)) +
  scale_color_manual(values = c("baetidae" = "black",
                                "heptageniidae" = "red"
  )) +
  labs(title = "Ll6") +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title =element_text(size = 12),
        legend.text = element_text(size = 10.5),
        axis.line = element_line(colour = "black", linewidth = 0.6))

sum.biplot

######
#Ll7

Ll7 <-subset(all3, community == "Ll7")
Ll7
Ll72 <- subset(Ll7, group == "heptageniidae" | group == "baetidae")
Ll72
#siber biplot Ll6
object <-createSiberObject(Ll72)
#setting up the plot
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot siber object
plotSiberObject(object,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'permille'),
                ylab = expression({delta}^15*N~'permille'),
                x.limits = c(-32, -21),
                y.limits = c(0, 5))

object$group.names
a.name.order <- c("Baetidae", "Heptageniidae")
legend("topright", a.name.order, 
       pch = c(1,1), col = c(1,2), lty = 1)

#Ll7 centroid plot
data.sum <- ddply(Ll72, c("group"), summarise,
                  iso1mn=mean(iso1),
                  iso1sd=sd(iso1),
                  iso2mn=mean(iso2),
                  iso2sd=sd(iso2))

data.sum
citation("plyr")
citation("ggplot2")

xlims <- aes(xmax = iso2mn + iso2sd, xmin = iso2mn - iso2sd)
ylims <- aes(ymax = iso1mn + iso1sd, ymin = iso1mn - iso1sd)

sum.biplot <- ggplot(data.sum, aes(x=iso2mn, y=iso1mn, colour = group)) +
  geom_point(size=3)+
  geom_errorbar(ylims, width=0.2)+
  geom_errorbarh(xlims, height=0.2) +
  ylab(expression(delta^{15}~N)) +
  xlab(expression(delta^{13}~C)) +
  scale_color_manual(values = c("baetidae" = "black",
                                "heptageniidae" = "red"
  )) +
  labs(title = "Ll7") +
  theme_light() +
  theme(axis.text = element_text(size = 12),
        axis.title =element_text(size = 12),
        legend.text = element_text(size = 10.5),
        axis.line = element_line(colour = "black", linewidth = 0.6))
sum.biplot



theme(text = element_text(size = 20),
      axis.text.x = element_text(angle = 90, hjust = 1)) 





#so we've got the plots, now we want to do the new overlaps and centroid distances
#Ll3
#lets start with overlap
#overlap
sea.overlap <- maxLikOverlap("1.leuctridae", "1.polycentropodidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

#leuctra and chloroperlidae
sea.overlap <- maxLikOverlap("Ll3.chloroperlidae", "Ll3.leuctridae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

#leuctra and perlodidae
sea.overlap <- maxLikOverlap("Ll3.chloroperlidae", "Ll3.limnephilidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

#leuctra nad limnephilidae

sea.overlap <- maxLikOverlap("Ll3.chloroperlidae", "Ll3.perlodidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

#leuctridae and tipulid

sea.overlap <- maxLikOverlap("Ll3.chloroperlidae", "Ll3.polycentropodidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)


#Polycentropod and chloroperlid

sea.overlap <- maxLikOverlap("Ll3.leuctridae", "Ll3.limnephilidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

#poly and perlodidae
sea.overlap <- maxLikOverlap("Ll3.leuctridae", "Ll3.perlodidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)


sea.overlap <- maxLikOverlap("Ll3.leuctridae", "Ll3.polycentropodidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)


sea.overlap <- maxLikOverlap("Ll3.limnephilidae", "Ll3.perlodidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)


sea.overlap <- maxLikOverlap("Ll3.limnephilidae", "Ll3.polycentropodidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

sea.overlap <- maxLikOverlap("Ll3.perlodidae", "Ll3.polycentropodidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)
library(geometry)
library(cxhull)
library(dplyr)
library(ggplot2)
library(SIBER)
library(rjags)

#now work out centroid distances for Ll3
#distance between centroids
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

ellipses.posterior <- siberMVN(object, parms, priors)
centroids <- siberCentroids(ellipses.posterior)
centroids
distances <- allCentroidVectors(centroids, upper = TRUE, do.plot = TRUE)
distances

#subset

Ll3.polycentropodidae.Ll3.chloroperlidae <- subset(distances, comparison == "Ll3.polycentropodidae.Ll3.chloroperlidae")
Ll3.polycentropodidae.Ll3.perlodidae <- subset(distances, comparison == "Ll3.polycentropodidae.Ll3.perlodidae")
Ll3.polycentropodidae.Ll3.limnephilidae <- subset(distances, comparison == "Ll3.polycentropodidae.Ll3.limnephilidae")
Ll3.polycentropodidae.Ll3.leuctridae <- subset(distances, comparison == "Ll3.polycentropodidae.Ll3.leuctridae")


Ll3.chloroperlidae.Ll3.perlodidae <- subset(distances, comparison == "Ll3.chloroperlidae.Ll3.perlodidae")
Ll3.chloroperlidae.Ll3.limnephilidae <- subset(distances, comparison == "Ll3.chloroperlidae.Ll3.limnephilidae")
Ll3.chloroperlidae.Ll3.leuctridae <-subset(distances, comparison == "Ll3.chloroperlidae.Ll3.leuctridae")

Ll3.perlodidae.Ll3.limnephilidae <- subset(distances, comparison == "Ll3.perlodidae.Ll3.limnephilidae")
Ll3.perlodidae.Ll3.leuctridae <- subset(distances, comparison == "Ll3.perlodidae.Ll3.leuctridae")

Ll3.limnephilidae.Ll3.leuctridae <- subset(distances, comparison == "Ll3.limnephilidae.Ll3.leuctridae")

summary(Ll3.polycentropodidae.Ll3.chloroperlidae)
summary(Ll3.polycentropodidae.Ll3.perlodidae)
summary(Ll3.polycentropodidae.Ll3.limnephilidae)
summary(Ll3.polycentropodidae.Ll3.leuctridae)

summary(Ll3.chloroperlidae.Ll3.perlodidae)
summary(Ll3.chloroperlidae.Ll3.limnephilidae)
summary(Ll3.chloroperlidae.Ll3.leuctridae)

summary(Ll3.perlodidae.Ll3.limnephilidae)
summary(Ll3.perlodidae.Ll3.leuctridae)

summary(Ll3.limnephilidae.Ll3.leuctridae)


#######
#Ll8 overlap
sea.overlap <- maxLikOverlap("Ll8.leuctridae", "Ll8.protonemura",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)
#leuctridae tipulid
sea.overlap <- maxLikOverlap("Ll8.leuctridae", "Ll8.tipulid",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)
#protonemura tipulid
sea.overlap <- maxLikOverlap("Ll8.protonemura", "Ll8.tipulid",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

#distance between centroids
#distance between centroids
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

ellipses.posterior <- siberMVN(object, parms, priors)
centroids <- siberCentroids(ellipses.posterior)
centroids
distances <- allCentroidVectors(centroids, upper = TRUE, do.plot = TRUE)
distances

#subset
Ll8.leuctridae.Ll8.protonemura <- subset(distances, comparison == "Ll8.leuctridae.Ll8.protonemura")
Ll8.leuctridae.Ll8.tipulid <- subset(distances, comparison == "Ll8.leuctridae.Ll8.tipulid")
Ll8.protonemura.Ll8.tipulid <- subset(distances, comparison =="Ll8.protonemura.Ll8.tipulid")

summary(Ll8.leuctridae.Ll8.protonemura)
summary(Ll8.leuctridae.Ll8.tipulid)
summary(Ll8.protonemura.Ll8.tipulid)

citation("rjags")



#####
#Ll6
sea.overlap <- maxLikOverlap("Ll6.heptageniidae", "Ll6.baetidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

#distance between centroids
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

ellipses.posterior <- siberMVN(object, parms, priors)
centroids <- siberCentroids(ellipses.posterior)
centroids
distances <- allCentroidVectors(centroids, upper = TRUE, do.plot = TRUE)
distances

#subset
Ll6.heptageniidae.Ll6.baetidae <- subset(distances, comparison =="Ll6.heptageniidae.Ll6.baetidae")

summary(Ll6.heptageniidae.Ll6.baetidae)



########
#Ll7
sea.overlap <- maxLikOverlap("Ll7.heptageniidae", "Ll7.baetidae",
                             object,
                             p.interval = 0.95,
                             n = 100,
                             do.plot = FALSE)

prop1 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.1"])
print(prop1)

prop2 <- as.numeric(sea.overlap["overlap"] / sea.overlap["area.2"])
print(prop2)

#distance between centroids
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

ellipses.posterior <- siberMVN(object, parms, priors)
centroids <- siberCentroids(ellipses.posterior)
centroids
distances <- allCentroidVectors(centroids, upper = TRUE, do.plot = TRUE)
distances

#subset
Ll7.heptageniidae.Ll7.baetidae <- subset(distances, comparison =="Ll7.heptageniidae.Ll7.baetidae")

summary(Ll7.heptageniidae.Ll7.baetidae)                         

