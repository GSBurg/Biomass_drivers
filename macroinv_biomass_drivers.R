# R Scripts for: Burgazzi et al. 2023 (the manuscript will be submitted soon)
# Title: Basin-scale variables drive macroinvertebrate biomass in low-order streams across different mountain ecoregions

# Authors: Gemma Burgazzi, Alex Laini, Pierluigi Viaroli, Stefano Fenoglio, Verena C. Schreiner, Ralf B. Schäfer, Alberto Doretto


# The code has been written by: 
# Dr. Gemma Burgazzi 
# Institute for Environmental Sciences, RPTU Kaiserslautern-Landau 
# Fortstraße 7, 76829 Landau in der Pfalz, Germany
# Email: gemmaburgazzi@gmail.com

# Revised by Dr. Verena Schreiner and Dr. Alex Laini

# R version 4.1.2 (2021-11-01) -- "Bird Hippie"

# With this study, we aim to:
# 1) evaluate the relationship between the biomass of macroinvertebrate communities (as mean biomass per organism per square meter) and several environmental variables
  # acting at three different spatial scales (patch, reach, and basin)
# 2) investigate transferability of the results among two geographic areas (the Maritime-Cottian Alps and the Tuscan-Emilian Apennine, hereafter ecoregions) by 
  # assessing the concordance of dominant drivers. 

# We sampled 10 streams in each geographic region, with 6 sampling point per stream

# Code sections:
# 1. Load libraries
# 2. Data import and manipulation
# 3. Create new predictors (with sparse PCA)
# 4. Mixed models
# 5. Variance partitioning
# 6. Figures

set.seed(667) #set the seed for reproducibility

#=============================== LOAD LIBRARIES
library(dplyr) # version 1.1.0, data manipulation
library(tibble) # version 3.1.8, data manipulation
library(pcaPP) # version 2.0.3, sparse PCA
library(lme4) # version 1.1.31, mixed modelling
library(partR2) # version 0.9.1, variance partitioning
library(ggplot2) # version 3.4.0, prepare figures
library(cowplot) # version 1.1.1, arrange figures

#library(vegan)
#source("multiplot.R")
#library(usdm) # version 1.1.18, check collinearity
#library(lmerTest) # version 3.1.3, mixed modelling

#=============================== DATA IMPORT and MANIPULATION

env <- read.csv("final_data_biomass.csv", h=T) # import the dataset

env$m_bio <- log((env$biomass/env$abu)/0.05) # compute mean macroinvertebrate per square meter
env <- env %>% mutate_if(is.character, as.factor) # convert all character variables to factor

env$stream <- factor(env$stream, levels=unique(env$stream)) # order streams alphabetically
env$sub_dom <- factor(env$sub_dom, levels = c("sand", "gravel", "cobbles", "boulders", "organic")) # order substrate categories according to size
str(env) # check the database

env_alp <- subset(env, ecoregion=="Alps") # subset the Alps database
env_ape <- subset(env, ecoregion=="Apennines") # subset the Apennines database
env_alp$stream <- droplevels(env_alp$stream) # remove extra levels
env_ape$stream <- droplevels(env_ape$stream) # remove extra levels


#=============================== CREATE NEW PREDICTORS (FROM sparse PCA)
# In the original dataset we have a high number of explanatory variables
# Moreover, the number of variables is not the same for each spatial scale
# To reduce dimensionality we use sparse PCA for env variables (one for each spatial scale and ecoregion)
# We extract PC1, PC2 and PC3 from each scale and build a new dataset of orthogonal predictors


## Sparse PCA for Patch scale

# Both for Alps and Apennines extract:
# Flow velocity (vel)
# Water depth (depth)
# Froude number (froude)
# Dominant substrate (sub_dom)
# Coarse Particulate Organic Matter (cpom)

# Alp
env_st_patch_alp <- subset(env_alp, select=c(vel, depth, froude, sub_dom, cpom)) # select patch scale variables
env_st_patch_alp <- model.matrix(~., env_st_patch_alp)[, -c(1,8)] #convert substrate into dummy variables. Also remove column 1 (intercept) and 8 (substrate organic, empty)
k.max <- 3 # set the number of wished axes
oTPO_D1_patch_alp <- opt.TPO(scale(env_st_patch_alp), k.max = k.max, method = "sd") #compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_patch_alp$pc)
oTPO_D1_patch_alp$pc$load
oTPO_D1_patch_alp$pc.noord$lambda
spc_patch_alp <- sPCAgrid(scale(env_st_patch_alp), k = k.max, lambda = oTPO_D1_patch_alp$pc.noord$lambda, method = "sd")
spc_patch_alp$loadings
spc_patch_alp_scores <- spc_patch_alp$scores%>%
  as.data.frame() %>%
  rownames_to_column("ID") # store SPC axes in the new db

#Ape
env_st_patch_ape <- subset(env_ape, select=c(vel, depth, froude, sub_dom, cpom)) # select patch scale variables
env_st_patch_ape <- model.matrix(~., env_st_patch_ape)[, -c(1)] # convert substrate into dummy variables. Also remove column 1 (intercept)
summary(env_st_patch_ape)
k.max <- 3 # set the number of wished axes
oTPO_D1_patch_ape <- opt.TPO(scale(env_st_patch_ape), k.max = k.max, method = "sd") # compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_patch_ape$pc)
oTPO_D1_patch_ape$pc$load
oTPO_D1_patch_ape$pc.noord$lambda
spc_patch_ape <- sPCAgrid(scale(env_st_patch_ape), k = k.max, lambda = oTPO_D1_patch_ape$pc.noord$lambda, method = "sd")
spc_patch_ape$loadings
spc_patch_ape_scores <- spc_patch_ape$scores%>%
  as.data.frame() %>%
  rownames_to_column("ID") # store SPC axes in the new db


## Sparse PCA for Reach scale

# Both for Alps and Apennines extract:
# Wetted channel width (wcw)
# Active channel width (acw)
# Water temperature (temp)
# pH (ph)
# Conductivity (cond_corr)
# Elevation above sea level (elevation)
# Strahler order (strahler)
# % of urban land use in 1 km buffer (Urb_1km)
# % of agricultural land use in 1 km buffer (Agr_1km)
# % of natural land use in 1 km buffer (Nat_1km)

#Alp
env_st_reach_alp <- subset(env_alp, select=c(wcw, acw, temp, ph, cond_corr, elevation, strahler, Urb_1km, Agr_1km, Nat_1km)) # select reach scale variables
k.max <- 3 # set the number of wished axes
oTPO_D1_reach_alp <- opt.TPO(scale(env_st_reach_alp), k.max = k.max, method = "sd") # compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_reach_alp$pc)
oTPO_D1_reach_alp$pc$load
oTPO_D1_reach_alp$pc.noord$lambda
spc_reach_alp <- sPCAgrid(scale(env_st_reach_alp), k = k.max, lambda = oTPO_D1_reach_alp$pc.noord$lambda, method = "sd")
spc_reach_alp$loadings
spc_reach_alp_scores <- spc_reach_alp$scores%>%
  as.data.frame() %>%
  rownames_to_column("ID") # store SPC axes in the new db

#Ape
env_st_reach_ape <- subset(env_ape, select=c(wcw, acw, temp, ph, cond_corr, elevation, strahler, Urb_1km, Agr_1km, Nat_1km)) # select reach scale variables
k.max <- 3 # set the number of wished axes
oTPO_D1_reach_ape <- opt.TPO(scale(env_st_reach_ape), k.max = k.max, method = "sd") # compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_reach_ape$pc)
oTPO_D1_reach_ape$pc$load
oTPO_D1_reach_ape$pc.noord$lambda
spc_reach_ape <- sPCAgrid(scale(env_st_reach_ape), k = k.max, lambda = oTPO_D1_reach_ape$pc.noord$lambda, method = "sd")
spc_reach_ape$loadings
spc_reach_ape_scores <- spc_reach_ape$scores%>%
  as.data.frame() %>%
  rownames_to_column("ID") # store SPC axes in the new db


## Sparse PCA for Basin scale

# Both for Alps and Apennines extract:
# Upslope basin area  (upslope_area)
# % of urban land use in the upslope basin (Urb)
# % of agricultural land use in the upslope basin (Agr)
# % of natural land use in the upslope basin (Nat)
# Mean value of cumulated daily precipitations (Prec_mean)
# Mean daily air temperature (Air_temp)

#Alp
env_st_basin_alp <- subset(env_alp, select=c(upslope_area, Urb, Agr, Nat, Prec_mean, Air_temp)) # select basin scale variables
k.max <- 3 # set the number of wished axes
oTPO_D1_basin_alp <- opt.TPO(scale(env_st_basin_alp), k.max = k.max, method = "sd") # compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_basin_alp$pc)
oTPO_D1_basin_alp$pc$load
oTPO_D1_basin_alp$pc.noord$lambda
spc_basin_alp <- sPCAgrid(scale(env_st_basin_alp), k = k.max, lambda = oTPO_D1_basin_alp$pc.noord$lambda, method = "sd")
spc_basin_alp$loadings
spc_basin_alp_scores <- spc_basin_alp$scores%>%
  as.data.frame() %>%
  rownames_to_column("ID") # store SPC axes in the new db

#Ape
env_st_basin_ape <- subset(env_ape, select=c(upslope_area, Urb, Agr, Nat, Prec_mean, Air_temp)) # select basin scale variables
k.max <- 3 # set the number of wished axes
oTPO_D1_basin_ape <- opt.TPO(scale(env_st_basin_ape), k.max = k.max, method = "sd") # compute a suggestion for lambda parameter (to be used in sPCAgrid)
summary(oTPO_D1_basin_ape$pc)
oTPO_D1_basin_ape$pc$load
oTPO_D1_basin_ape$pc.noord$lambda
spc_basin_ape <- sPCAgrid(scale(env_st_basin_ape), k = k.max, lambda = oTPO_D1_basin_ape$pc.noord$lambda, method = "sd")
spc_basin_ape$loadings
spc_basin_ape_scores <- spc_basin_ape$scores %>%
  as.data.frame() %>%
  rownames_to_column("ID") # store SPC axes in the new db

# rename SPC variables
colnames(spc_patch_alp_scores)[2:4] <- colnames(spc_patch_ape_scores)[2:4] <- paste(c("PC1", "PC2", "PC3"), "patch", sep = "_")
colnames(spc_reach_alp_scores)[2:4] <- colnames(spc_reach_ape_scores)[2:4] <- paste(c("PC1", "PC2", "PC3"), "reach", sep = "_")
colnames(spc_basin_alp_scores)[2:4] <- colnames(spc_basin_ape_scores)[2:4] <- paste(c("PC1", "PC2", "PC3"), "basin", sep = "_")

# add the new SPC variables to the Alps database
env_alp <- env_alp %>% 
  rownames_to_column("ID") %>%
  inner_join(.,spc_patch_alp_scores)%>%
  inner_join(.,spc_reach_alp_scores)%>%
  inner_join(.,spc_basin_alp_scores)

# add the new SPC variables to the Apennines database
env_ape <- env_ape %>% 
  rownames_to_column("ID") %>%
  inner_join(.,spc_patch_ape_scores)%>%
  inner_join(.,spc_reach_ape_scores)%>%
  inner_join(.,spc_basin_ape_scores)



#=============================== MIXED MODELS with SPC axes
# Both for Alps and Apennines we include 3 SPC for each spatial scale as fixed factor, whereas the stream is included as random factor


## Alps
alp.lme <- lmer(m_bio ~ PC1_patch+PC2_patch+PC3_patch +
                  PC1_reach+PC2_reach+PC3_reach +
                  PC1_basin+PC2_basin+PC3_basin +
                  (1|stream), data = env_alp) 

anova(alp.lme, type = "II") # summarize the model with Type II anova

# Model check
plot(alp.lme)
text(qqnorm(resid(alp.lme)))
qqline(resid(alp.lme))
hist(resid(alp.lme))


## Apennines
ape.lme <- lmer(m_bio ~ PC1_patch+PC2_patch+PC3_patch +
                  PC1_reach+PC2_reach+PC3_reach +
                  PC1_basin+PC2_basin+PC3_basin +
                  (1|stream), data = env_ape)

anova(ape.lme, type = "II") # summarize the model with Type II anova

# Model check
plot(ape.lme)
text(qqnorm(resid(ape.lme)))
qqline(resid(ape.lme))
hist(resid(ape.lme))



#=============================== VARIANCE PARTITIONING on MIXED MODELS
# Run variance partitioning on the two mixed models previously created
# Be careful, partR2 is time consuming. For a quick try-out reduce the nboot

# Alps
alp_vp <- partR2(alp.lme, partbatch = list(Patch = c("PC1_patch","PC2_patch","PC3_patch"), 
                                           Reach = c("PC1_reach","PC2_reach","PC3_reach"),
                                           Basin = c("PC1_basin","PC2_basin","PC3_basin")), nboot=1000)
summary(alp_vp)

# Apennines
ape_vp <- partR2(ape.lme, partbatch = list(Patch = c("PC1_patch","PC2_patch","PC3_patch"), 
                                           Reach = c("PC1_reach","PC2_reach","PC3_reach"),
                                           Basin = c("PC1_basin","PC2_basin","PC3_basin")), nboot=1000)
summary(ape_vp)


#=============================== MANUSCRIPT FIGURES

## Figure 1
# Figure 1 was created with Quantum GIS

## Figure 2
g1 <- ggplot(env, aes(x=ecoregion, y=m_bio, fill=ecoregion)) + geom_boxplot(color="black", linewidth=0.8) + 
  scale_fill_manual(values = alpha(c("blue", "darkgreen"), 0.7)) + 
  theme_test()+ 
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold", colour = "black"), axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12,face="bold")) + 
  ggtitle("(a) Biomass variability among ecoregions") + ylab("log-transformed biomass/abundance ratio") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("") 
g1


g2 <- ggplot(env, aes(y=stream, x=m_bio, fill=ecoregion)) + geom_boxplot(color="black", linewidth=0.8) + 
  scale_fill_manual(values = alpha(c("blue", "darkgreen"), 0.7)) + 
  theme_test()+ 
  theme(legend.position = "none", plot.title = element_text(size = 14, face = "bold", colour = "black"), axis.text=element_text(size=12, colour = "black"), 
        axis.title=element_text(size=12,face="bold")) + 
  ggtitle("(b) Biomass variability among streams") + ylab("log-transformed biomass/abundance ratio") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("") 
g2


tiff("Figure_2.tiff", width = 12, height = 6, units = 'in', res = 600, compression = 'none')
#par(mfrow=c(1,2), mai=c(0.3, 0.3, 1, 0.3))
plot_grid(g1, g2, ncol = 2, nrow = 1)
dev.off()


## Figure 3

# plots of the varparts (one for each ecoregion)
p_alp <- forestplot(alp_vp, type = "R2") + theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11)) +
  ggtitle("(a) Alps")

p_ape <- forestplot(ape_vp, type = "R2")  + theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 11)) +
  ggtitle("(b) Apennines")

tiff("Figure_3.tiff", width = 10, height = 4.5, units = 'in', res = 600, compression = 'none')
par(mfrow=c(1,2), mai=c(0.3, 0.3, 1, 0.3))
plot_grid(p_alp, p_ape, ncol = 2, nrow = 1)
dev.off()



