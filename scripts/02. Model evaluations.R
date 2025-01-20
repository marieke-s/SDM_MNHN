<<<<<<< HEAD
################################################################################
# Species distribution modelling with biomod
#
# AUTHOR: Boris Leroy
# 
# LICENSE: GPL v3
# GPL-3.0-or-later  
#
#  - Permissions 	 	
# 
# Commercial use
# Distribution
# Modification
# Patent use
# Private use
# 
#  - Conditions
# Disclose source
# License and copyright notice
# Same license
# State changes
# 
#  - Limitations
# Liability
# Warranty
#
# Anyone can copy, modify and distribute this code. You have to credit the 
# author and include the license and copyright notice with each and every 
# distribution. Any modifications of this code base MUST be distributed with the
# same license, GPLv3.
#
################################################################################


library(biomod2)
library(ggplot2)
library(plyr)

sp_list <- read.csv("data_cours/species_list.csv", sep = ";")


evals_allsp <- data.frame()

for (i in 1:nrow(sp_list))
{
  # 1. On charge les modèles calibrés et leurs évals pour l'espèce actuelle
  sp <- sp_list$sp[i]
  cat(paste("----", Sys.time(), sp, " evaluation initialised ----\n", sep = " "))

  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))

  
  # 2. On récupère les évaluations ROC et TSS
  evals_allsp <- rbind.fill(evals_allsp,
                            data.frame(Species = sp,
                                       get_evaluations(model_runs)))
}


ggplot(evals_allsp, aes(y = validation,
                        x = Species,
                        col = algo)) +
  geom_boxplot() +
  facet_wrap(~metric.eval, scales = "free")
  
ggplot(evals_allsp, aes(y = validation,
                        x = Species)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(metric.eval ~ algo, scales = "free") +
  coord_flip()



levels(as.factor(evals_allsp$Species))


ggplot(evals_allsp) + 
  geom_boxplot(aes(x = algo, y = validation, col = Species)) + 
  facet_grid(metric.eval ~ ., scale = "free_y")



# Fonctions graphiques intégrées à biomod2 :

model_runs <- readRDS("models/Tyrion.lannisterii/model_runs.RDS")
bm_PlotEvalMean(model_runs)

=======
################################################################################
# Species distribution modelling with biomod
#
# AUTHOR: Boris Leroy
# 
# LICENSE: GPL v3
# GPL-3.0-or-later  
#
#  - Permissions 	 	
# 
# Commercial use
# Distribution
# Modification
# Patent use
# Private use
# 
#  - Conditions
# Disclose source
# License and copyright notice
# Same license
# State changes
# 
#  - Limitations
# Liability
# Warranty
#
# Anyone can copy, modify and distribute this code. You have to credit the 
# author and include the license and copyright notice with each and every 
# distribution. Any modifications of this code base MUST be distributed with the
# same license, GPLv3.
#
################################################################################


library(biomod2)
library(ggplot2)
library(plyr)

sp_list <- read.csv("data_cours/species_list.csv", sep = ";")


evals_allsp <- data.frame()

for (i in 1:nrow(sp_list))
{
  # 1. On charge les modèles calibrés et leurs évals pour l'espèce actuelle
  sp <- sp_list$sp[i]
  cat(paste("----", Sys.time(), sp, " evaluation initialised ----\n", sep = " "))

  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))

  
  # 2. On récupère les évaluations ROC et TSS
  evals_allsp <- rbind.fill(evals_allsp,
                            data.frame(Species = sp,
                                       get_evaluations(model_runs)))
}


ggplot(evals_allsp, aes(y = validation,
                        x = Species,
                        col = algo)) +
  geom_boxplot() +
  facet_wrap(~metric.eval, scales = "free")
  
ggplot(evals_allsp, aes(y = validation,
                        x = Species)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(metric.eval ~ algo, scales = "free") +
  coord_flip()



levels(as.factor(evals_allsp$Species))


ggplot(evals_allsp) + 
  geom_boxplot(aes(x = algo, y = validation, col = Species)) + 
  facet_grid(metric.eval ~ ., scale = "free_y")



# Fonctions graphiques intégrées à biomod2 :

model_runs <- readRDS("models/Tyrion.lannisterii/model_runs.RDS")
bm_PlotEvalMean(model_runs)

>>>>>>> 5768564 (2nd commit)
