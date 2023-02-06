#### SCRIPT PAPER BIOMASS PRIN (final version version) ####
# Contributors: G. Burgazzi, R.B. Schaefer, A. Doretto, V. Schreiner
set.seed(667)

#=============================== LOAD LIBRARIES
library(ggplot2)
library(cowplot)
library(lme4)
library(lmerTest)
library(vegan)
library(biomonitoR)
library(usdm)
library(partR2)
library(pcaPP)
library(dplyr)
source("multiplot.R")



#=============================== DATA IMPORT and MANIPULATION
env <- read.table("biomass_env.txt", h=T)
taxa <- read.table("biomass_taxa.txt", h=T)
env$m_bio <- log((env$biomass/env$abu)/0.05)
env$stream <- factor(env$stream, levels=unique(env$stream))
env$sub_dom <- as.factor(env$sub_dom)
unique(env$sub_dom)
env$sub_dom <- factor(env$sub_dom, levels = c("sand", "gravel", "cobbles", "boulders", "organic"))
str(env)

env_alp <- subset(env, env$ecoregion=="Alps")
env_ape <- subset(env, env$ecoregion=="Apennines")

# The Alpine dataset has more than 6 samples per stream. We randomly select 6 from each stream to have a balanced dataset
env_alp <- env_alp %>% group_by(stream) %>% slice_sample(n = 6)
env <- rbind(env_alp, env_ape)


#=============================== CREATE NEW PREDICTORS (FROM sparse PCA)
# See: https://github.com/rbslandau/graf_spidertraits
#### sparse PCA for env variables (one for each spatial scale)
# extract PC1, PC2 and PC3 from each scale and build a new dataset of orthogonal predictors

## create an empty db for storing the new variables
new_cov <- data <- data.frame(matrix(NA, nrow = nrow(env), ncol = 9))
rownames(new_cov) <- rownames(env)
colnames(new_cov) <- c(paste(c("PC1", "PC2", "PC3"), "patch", sep = "_"),
                       paste(c("PC1", "PC2", "PC3"), "reach", sep = "_"),
                       paste(c("PC1", "PC2", "PC3"), "basin", sep = "_"))


## Sparse PCA for Patch scale
#Alp
env_st_patch_alp <- subset(env_alp, select=c(vel, depth, froude, sub_dom, cpom)) #select patch scale variables
env_st_patch_alp <- model.matrix(~., env_st_patch_alp) #convert substrate into dummy variables
summary(env_st_patch_alp)
env_st_patch_alp <- subset(env_alp, select=c(vel, depth, froude, sub_dom, cpom)) #select patch scale variables
env_st_patch_alp <- model.matrix(~., env_st_patch_alp)[, -c(1,8)] #convert substrate into dummy variables
summary(env_st_patch_alp)
k.max <- 3 # set the number of wished axes
oTPO_D1_patch_alp <- opt.TPO(scale(env_st_patch_alp), k.max = k.max, method = "sd") #compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_patch_alp$pc)
oTPO_D1_patch_alp$pc$load
oTPO_D1_patch_alp$pc.noord$lambda
spc_patch_alp <- sPCAgrid(scale(env_st_patch_alp), k = k.max, lambda = oTPO_D1_patch_alp$pc.noord$lambda, method = "sd")
spc_patch_alp$loadings
new_cov[1:60,1:3] <- spc_patch_alp$scores #store PC1:3 in the new db

#Ape
env_st_patch_ape <- subset(env_ape, select=c(vel, depth, froude, sub_dom, cpom)) #select patch scale variables
env_st_patch_ape <- model.matrix(~., env_st_patch_ape)[, -c(1)] #convert substrate into dummy variables
summary(env_st_patch_ape)
k.max <- 3 # set the number of wished axes
oTPO_D1_patch_ape <- opt.TPO(scale(env_st_patch_ape), k.max = k.max, method = "sd") #compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_patch_ape$pc)
oTPO_D1_patch_ape$pc$load
oTPO_D1_patch_ape$pc.noord$lambda
spc_patch_ape <- sPCAgrid(scale(env_st_patch_ape), k = k.max, lambda = oTPO_D1_patch_ape$pc.noord$lambda, method = "sd")
spc_patch_ape$loadings
new_cov[61:120,1:3] <- spc_patch_ape$scores #store PC1:3 in the new db


## Sparse PCA for Reach scale
#Alp
env_st_reach_alp <- subset(env_alp, select=c(wcw, acw, temp, ph, cond_corr, altitude, strahler, Urb_1km, Agr_1km, Nat_1km)) #select reach scale variables
k.max <- 3 # set the number of wished axes
oTPO_D1_reach_alp <- opt.TPO(scale(env_st_reach_alp), k.max = k.max, method = "sd") #compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_reach_alp$pc)
oTPO_D1_reach_alp$pc$load
oTPO_D1_reach_alp$pc.noord$lambda
spc_reach_alp <- sPCAgrid(scale(env_st_reach_alp), k = k.max, lambda = oTPO_D1_reach_alp$pc.noord$lambda, method = "sd")
spc_reach_alp$loadings
new_cov[1:60,4:6] <- spc_reach_alp$scores #store PC1:3 in the new db

#Ape
env_st_reach_ape <- subset(env_ape, select=c(wcw, acw, temp, ph, cond_corr, altitude, strahler, Urb_1km, Agr_1km, Nat_1km)) #select reach scale variables
k.max <- 3 # set the number of wished axes
oTPO_D1_reach_ape <- opt.TPO(scale(env_st_reach_ape), k.max = k.max, method = "sd") #compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_reach_ape$pc)
oTPO_D1_reach_ape$pc$load
oTPO_D1_reach_ape$pc.noord$lambda
spc_reach_ape <- sPCAgrid(scale(env_st_reach_ape), k = k.max, lambda = oTPO_D1_reach_ape$pc.noord$lambda, method = "sd")
spc_reach_ape$loadings
new_cov[61:120,4:6] <- spc_reach_ape$scores #store PC1:3 in the new db


# Sparse PCA for Basin scale
#Alp
env_st_basin_alp <- subset(env_alp, select=c(upslope_area, Urb, Agr, Nat, Prec_mean, Air_temp)) #select basin scale variables
k.max <- 3 # set the number of wished axes
oTPO_D1_basin_alp <- opt.TPO(scale(env_st_basin_alp), k.max = k.max, method = "sd") #compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_basin_alp$pc)
oTPO_D1_basin_alp$pc$load
oTPO_D1_basin_alp$pc.noord$lambda
spc_basin_alp <- sPCAgrid(scale(env_st_basin_alp), k = k.max, lambda = oTPO_D1_basin_alp$pc.noord$lambda, method = "sd")
spc_basin_alp$loadings
new_cov[1:60,7:9] <- spc_basin_alp$scores #store PC1:3 in the new db

#Ape
env_st_basin_ape <- subset(env_ape, select=c(upslope_area, Urb, Agr, Nat, Prec_mean, Air_temp)) #select basin scale variables
k.max <- 3 # set the number of wished axes
oTPO_D1_basin_ape <- opt.TPO(scale(env_st_basin_ape), k.max = k.max, method = "sd") #compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_basin_ape$pc)
oTPO_D1_basin_ape$pc$load
oTPO_D1_basin_ape$pc.noord$lambda
spc_basin_ape <- sPCAgrid(scale(env_st_basin_ape), k = k.max, lambda = oTPO_D1_basin_ape$pc.noord$lambda, method = "sd")
spc_basin_ape$loadings
new_cov[61:120,7:9] <- spc_basin_ape$scores #store PC1:3 in the new db



# add the new variables to the env db
env <- cbind(env, new_cov)
env_alp <- subset(env, env$ecoregion=="Alps")
env_ape <- subset(env, env$ecoregion=="Apennines")




#============================= VARPART on MIXED MODELS with PCs
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8162244/


## Alps

# mixed model
alp.lme <- lmer(m_bio ~ PC1_patch+PC2_patch+PC3_patch +
                  PC1_reach+PC2_reach+PC3_reach +
                  PC1_basin+PC2_basin+PC3_basin +
                  (1|stream), data = env_alp[-12,])


plot(alp.lme)
text(qqnorm(resid(alp.lme)))
qqline(resid(alp.lme))
anova(alp.lme, type = "II")

# varpart
alp_vp <- partR2(alp.lme, partbatch = list(Patch = c("PC1_patch","PC2_patch","PC3_patch"), 
                                           Reach = c("PC1_reach","PC2_reach","PC3_reach"),
                                           Basin = c("PC1_basin","PC2_basin","PC3_basin")), nboot=1000)
summary(alp_vp)



## Apennines

# mixed model
ape.lme <- lmer(m_bio ~ PC1_patch+PC2_patch+PC3_patch +
                  PC1_reach+PC2_reach+PC3_reach +
                  PC1_basin+PC2_basin+PC3_basin +
                  (1|stream), data = env_ape[-c(53),])

plot(ape.lme)
text(qqnorm(resid(ape.lme)))
qqline(resid(ape.lme))
anova(ape.lme, type = "II")

# varpart
ape_vp <- partR2(ape.lme, partbatch = list(Patch = c("PC1_patch","PC2_patch","PC3_patch"), 
                                           Reach = c("PC1_reach","PC2_reach","PC3_reach"),
                                           Basin = c("PC1_basin","PC2_basin","PC3_basin")), nboot=1000)
summary(ape_vp)



# plots of the varparts
p_alp <- forestplot(alp_vp, type = "R2") + theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11)) +
  ggtitle("a. Alps")

p_ape <- forestplot(ape_vp, type = "R2")  + theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11)) +
  ggtitle("b. Apennines")


multiplot(p_alp, p_ape, cols = 2)




#=============================== Manuscript figures

## Figure 2
g1 <- ggplot(env, aes(x=ecoregion, y=m_bio, fill=ecoregion)) + geom_boxplot() + 
  scale_fill_manual(values = c("blue", "darkgreen")) + 
  theme_test()+ 
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold")) + 
  ggtitle("a. Biomass") + ylab("log-transformed biomass/abundance ratio") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("")
g1


g2 <- ggplot(env, aes(y=stream, x=m_bio, fill=ecoregion)) + 
  stat_summary(geom="bar", fun=mean) +
  stat_summary(geom="errorbar", fun.data = mean_se, size=0.8, width = 0.25) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(size = 8, colour = "black", face= "bold", angle = 0, vjust = 0.5, hjust=1), axis.title.y=element_blank()) + 
  ggtitle("b. Biomass variability among streams") + xlab("log-transformed biomass/abundance ratio") + 
  theme(plot.title = element_text(size = 14, face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = c("blue", "darkgreen"))
g2 #not included in the manuscript


### Size check
# Prepare data for the analysis.
taxanames <- colnames(taxa)
macro <- as.data.frame(t(taxa))
macro$Taxa <- taxanames
data_bio <- as_biomonitor(macro, traceB = TRUE)
data_agr <- aggregate_taxa(data_bio)

# Functional metrics
data_ts <- assign_traits(data_agr)
data_ts_av <- average_traits( data_ts )
col_blocks <- c( 8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8 )


#CWM (community weighted matrix) for the maximum size
tts <- t(taxa)
tts <- tts[!(row.names(tts) %in% c("Mermithidae", "Irudinea", "Niphargidae", "Tricladida","Trombidiformes")),] #I had to remove these taxa cause there are no traits for them in biomonitoR dataset

size <- paste("SIZE", 1:7, sep ="_")
size <- data_ts_av[,size]
rownames(size) <- data_ts_av$Taxa
colnames(size) <- c("<0.25 cm","0.25-0.5 cm",	"0.5-1 cm",	"1-2 cm",	"2-4 cm",	"4-8 cm",	">8 cm")

temp <- as.data.frame(matrix(data = NA, nrow = nrow(tts), ncol = 7))
rownames(temp) <- rownames(tts)
colnames(temp) <- colnames(size)
temp2 <- as.data.frame(matrix(data = NA, nrow = nrow(env), ncol = 7))
rownames(temp2) <- rownames(env)
colnames(temp2) <- colnames(size)

for (j in 1:nrow(env)) {
  for (i in 1:nrow(tts)) {
    temp[i, ] <-  tts[i,j]*size[i,]
  }
  temp2[j,] <- apply(temp, 2, sum)
}

cwm <- temp2
cwm$Ecoregion <- env$ecoregion
library(reshape2)
cwm_long <- melt(cwm, id.vars = c("Ecoregion"), variable.name = "size")

b1 <- ggplot(cwm_long, aes(x=size, y=log(value + 1), fill=Ecoregion)) + geom_boxplot() +
  theme_test() +
  scale_fill_manual(values = c("blue", "darkgreen")) +
  ylab("log-transformed CWM for size classes") +
  xlab("Size classes") +
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("b. Size classes") 
b1



tiff("Fig_2.tiff", width = 10, height = 5.5, units = 'in', res = 600, compression = 'none')
ggdraw() +
  draw_plot(g1, x = 0, y = 0, width = 0.35, height = 1) +
  draw_plot(b1, x = 0.4, y = 0, width = 0.6, height = 1)
dev.off()

## Figure 3
tiff("Fig_3_new.tiff", width = 10, height = 4.5, units = 'in', res = 600, compression = 'none')
par(mfrow=c(1,2), mai=c(0.3, 0.3, 1, 0.3))
multiplot(p_alp, p_ape, cols = 2)
dev.off()


