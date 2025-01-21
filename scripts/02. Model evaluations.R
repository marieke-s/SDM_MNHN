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
# Load libraries
library(biomod2)
library(ggplot2)
library(plyr)

# Load species list
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")

# Df pour récupérer les évaluations
evals_allsp <- data.frame()

# Récupération des évaluations
for (i in 1:nrow(sp_list))
{
  # Test loop
  # i <- 1 
  
  # 1. On charge les modèles calibrés et leurs évals pour l'espèce actuelle
  sp <- sp_list$sp[i]
  cat(paste("----", Sys.time(), sp, " evaluation initialised ----\n", sep = " "))

  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))

  
  # 2. On récupère les évaluations ROC et TSS et on y ajoute l'espèce --> ce df va permettre de faire du ggplot après
  evals_allsp <- rbind.fill(evals_allsp,
                            data.frame(Species = sp,
                                       get_evaluations(model_runs)))
  # get_evaluations(model_runs) : cutoff : à ignorer car Boyce n'a pas de seuil normalement
  # 
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

# Interpretation Boyce, ROC, TSS
# AUC = ROC : calcul qui ressemble bcp au TSS et indice de Jaccard
# Plus les spp sont rares plus les indicateurs sont élévés
# ROC > 0.8 = ok > 0.9 = excellent
# Boyce, ROC et TSS n PA avec un gros déséquilibre : pas fiables (<10% de présences par rapport aux absences) --> il faut utiliser les l'indice de Jaccard

# Probème avec TSS et ROC : basés sur une matrice de confusion, prennent en compte les vrais négatifs, qui sont artificielllement gonglés quand on a beaucoup d'absences --> on a bcp de TN et donc un TSS élevé : les faux positifs ne sont donc pas pénalisés. --> d'où biais mathématiques pour les espèces avec une petite aire de répartition / PA déséquilibrée.

# Solution : utiliser l'indice de Jaccard 
# On évalue le chevauchement entre l'aire observée vs l'aire prédite (cf schema carnet de note rouge page 48)






