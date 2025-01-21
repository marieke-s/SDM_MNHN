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
library(terra)

if(!dir.exists("outputs"))
{
  dir.create("outputs", recursive = T)
}
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")


proj_names <- readRDS("data/projection_names.RDS")

stacks <- c("baseline",
            apply(proj_names[, c("gcm_lu", "scenar_clim", "year_lu")], 1, 
                  paste, collapse = "_"))


scenarios <- c("ssp245_2070", "ssp585_2070")

allsp_cutoffs <- readRDS("data/jaccard_cutoffs.RDS")



#### NOTE ####
# This script creates a lot of ensemble projections, but YOU DO NOT NEED TO
# USE ALL OF THEM. Choose only the ones that you need for your specific study.
# Conceptualise your study BEFORE using this script and decide what you want
# to use.


# 7 stacks will be created:  

# - em_cont  # Ensemble model suitabilities
# - em_sd # Ensemble model uncertainty (sd)
# - em_bin  # Ensemble model potential P/A
# - SRC_stack # Species Range Change stack based on em_bin

# - em_penalised_cont # Penalised EM suitabilities
# - em_penalised_bin # Penalised EM potential P/A


# - em_CA # EM Committee Averaging



# Note: a penalised SRC stack is not computed, as it does not make sense:
# Current values are based on a single climate projection
# Future values are based on 3 climate projection
# Hence, the standard deviation in the future is increased by the use of different GCM
# Comparing current values to future values would overestimate the change because of 
# the increased variability between models + GCMs


for (i in 1:nrow(sp_list))
{
  sp <- sp_list$sp[i]
  cat(paste("----", Sys.time(), sp, "stack creation initialised ----\n", sep = " "))
  
  
  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))
  # Récupération des données entrées (notamment pour avoir les coordonnées 
  # des pts de présence)
  input_data <- get_formal_data(model_runs)
  sp_coords <- input_data@coord
  
  
  # Récupération des seuils issus de Jaccard
  tmp <- reshape2::melt(allsp_cutoffs[sp, , , ])
  cur_cutoffs <- tmp$value
  # On va s'assurer que les noms pour les cutoffs correspondent bien aux noms
  # des projections faites par biomod
  names(cur_cutoffs) <- paste(sp, tmp$pa.run, tmp$cv.run, tmp$model, sep = "_")
  
  ##### Partie 1 : baseline #####
  # Cartes de probabilité de présence
  cur_stack <- rast(paste0("models/", sp, "/proj_", stacks[1],
                           "/proj_", stacks[1], "_", sp, ".tif"))
  # Transformation en présence-absence par le seuil de Jaccard
  cur_binary_stack <- bm_BinaryTransformation(cur_stack,
                                              cur_cutoffs[names(cur_stack)])
  
  # Ensemble model: probabilités
  cur_em <- mean(cur_stack) # Ensemble model
  cur_sd <- app(cur_stack, sd) # Uncertainty (standard deviation)
  cur_em_penalised <- cur_em - cur_sd # Penalised probability of presence
  cur_em_penalised[cur_em_penalised < 0] <- 0
  
  # Ensemble model: conversion en présence-absence
  # Ici j'illustre avec l'indice de Jaccard
  # Mais il est facile de reprendre plutôt le TSS ou ROC 
  # en allant chercher les évaluations dans l'objet d'EM de biomod
  
  
  obs_data <- input_data@data.species
  obs_data[is.na(obs_data)] <- 0 # Pour les présences pseudoabs,
  # On remplace les pseudoabs par des absences
  # Notez bien que c'est une hypothèse très forte (très fausse), mais 
  # c'est la seule possibilité en presence pseudoabs
  # A éviter si vous pouvez vous en passer. Soyez inventifs sur votre 
  # utilisation des sorties continues des modèles !
  # Cf. e.g. la méthode que j'ai utilisé dans le rapport SDMs Corse
  # https://borisleroy.com/sdms-pna-corse/
  
  # On extrait les probas aux coordonnées de présence
  # Notez qu'on le fait sur tous les points de présence (contrairement au script
  # Jaccard où on ne le faisait que sur les données d'évaluation)
  # Ce n'est pas un problème ici car **on n'évalue pas le modèle** : 
  # on recherche simplement le seuil optimal
  pred_data <- extract(cur_em, sp_coords)
  jaccard_test <- NULL
  for(cutoff in seq(0, 1000, by = 1))
  {
    pred.pa <- pred_data$mean
    pred.pa[pred.pa < cutoff] <- 0
    pred.pa[pred.pa >= cutoff] <- 1
    TP <- length(which(obs_data == 1 & pred.pa == 1))
    FN <- length(which(obs_data == 1 & pred.pa == 0))
    FP <- length(which(obs_data == 0 & pred.pa == 1))
    jaccard <- TP / (TP + FP + FN)
    jaccard_test <- rbind.data.frame(jaccard_test,
                                     data.frame(cutoff = cutoff,
                                                TP = TP,
                                                FN = FN,
                                                FP = FP,
                                                jaccard = jaccard))
  }
  em_cutoff <- mean(jaccard_test$cutoff[
    which(jaccard_test$jaccard == max(jaccard_test$jaccard))
  ])
  
  # Ensemble model transformé en présence-absence
  cur_em_binary <- bm_BinaryTransformation(cur_em, em_cutoff)
  # Transformation en présence-absence des probas pénalisées par l'incertitude
  cur_em_binary_penalised <- bm_BinaryTransformation(cur_em_penalised, 
                                                     em_cutoff)
  # Committee averaging
  cur_em_CA <- sum(cur_binary_stack)
  cur_em_CA[cur_em_CA < (nlyr(cur_binary_stack) / 2)] <- 0
  cur_em_CA[cur_em_CA >= (nlyr(cur_binary_stack) / 2)] <- 1
  
  
  # On prépare nos stacks qui contiendront tous les rasters finaux ici
  em_cont <- cur_em # Ensemble model suitabilities
  em_sd <- cur_sd # Ensemble model uncertainty (sd)
  em_bin <- cur_em_binary # Ensemble model potential P/A
  em_penalised_cont <- cur_em_penalised # Penalised EM suitabilities
  em_penalised_bin <- cur_em_binary_penalised # Penalised EM potential P/A
  em_CA <- cur_em_CA # EM Committee Averaging
  
  
  ##### Partie 2 : scénarios futurs #####
  # La difficulté ici est qu'il faut réunir les projections des différents GCMs,
  # stockés par biomod dans différents objets 
  for(s in scenarios)
  {
    cat(paste("\n-- Calculations for scenario:", s, "\n", sep = " "))
    # Il faut tout d'abord réunir les résultats de tous les
    # GCMs de notre scénario dans un seul et unique stack
    subStacks <- stacks[grep(s, stacks)]
    for (s1 in subStacks)
    {
      if(s1 == subStacks[1])
      {
        future_stack <- rast(paste("models/", sp, "/proj_", s1, "/", 
                                   "proj_", s1, "_", sp, ".tif", sep =""))
        order_list <- names(future_stack) # On récupère les noms que biomod a
        # donné, afin de garder le lien entre modèles et cutoffs à appliquer
        
        # Ensuite on renomme les couches en ajoutant le nom du GCM au début
        names(future_stack) <-  paste(s1, names(future_stack), sep = "_")
      } else
      {
        tmp <- rast(paste("models/", sp, "/proj_", s1, "/",
                          "proj_", s1, "_", sp, ".tif", sep =""))
        names(tmp) <- paste(s1, names(tmp), sep = "_")
        future_stack <- c(future_stack,
                          tmp)
      }
    }
    ##### Attention: si votre protocole est trop lourd (trop de modèles, 
    # trop de répétitions, trop de GCMs, taille de la zone d'étude trop grande 
    # ou résolution trop fine), les étapes suivantes pourraient être
    # trop lourdes pour un pc moyen. Dans ce cas il faudra
    # se "contenter" de faire la moyenne entre les modèles 
    # d'ensemble projetés par biomod sur les différents GCM
    
    # Cartes de présence-absence pour chaque modèle du futur
    future_binary_stack <- bm_BinaryTransformation(future_stack, 
                                                   rep(cur_cutoffs[order_list], 
                                                       length(subStacks)))
    
    # Ensemble model: probabilités
    future_em <- mean(future_stack) # Ensemble model
    future_sd <- app(future_stack, sd) # Uncertainty (standard deviation)
    future_em_penalised <- future_em - future_sd # Penalised probability of presence
    future_em_penalised[future_em_penalised < 0] <- 0
    
    # Ensemble model transformé en présence-absence
    future_em_binary <- bm_BinaryTransformation(future_em, em_cutoff)
    
    # Transformation en présence-absence des probas pénalisées par l'incertitude
    future_em_binary_penalised <- bm_BinaryTransformation(future_em_penalised, 
                                                          em_cutoff)
    
    # Committee averaging
    future_em_CA <- sum(future_binary_stack)
    future_em_CA[future_em_CA < (nlyr(future_binary_stack) / 2)] <- 0
    future_em_CA[future_em_CA >= (nlyr(future_binary_stack) / 2)] <- 1
    
    # Range size change (basé sur les rasters NON PENALISES)
    # WARNING: ne pas utiliser les rasters pénalisés ici, car
    # l'incertitude est plus grande dans le futur à cause de 
    # l'utilisation de différents GCMs
    # De toute façon à ce stade j'espère vous avoir convaincu de favoriser 
    # la proba plutôt que les stacks binaires !!!
    SRCRaster <- future_em_binary - 2 * cur_em_binary
    # -2 : aire perdue
    # -1 : aire stable
    # 0 : aire inoccupée dans les 2 périodes
    # 1 : aire nouvelle
    
    # tmp <- SRCRaster
    # tmp[tmp != -2] <- 0
    # tmp[tmp == -2] <- 1
    # global(tmp * cellSize(tmp), fun = "sum", na.rm = TRUE)
    

    # On complète nos stacks avec les résults de l'EM pour le scénario actuel
    em_cont <- c(em_cont,
                 future_em) # Ensemble model suitabilities
    em_sd <- c(em_sd,
               future_sd) # Ensemble model uncertainty (sd)
    em_bin <- c(em_bin,
                future_em_binary) # Ensemble model potential P/A
    em_penalised_cont <- c(em_penalised_cont,
                           future_em_penalised) # Penalised EM suitabilities
    em_penalised_bin <- c(em_penalised_bin,
                          future_em_binary_penalised) # Penalised EM potential P/A
    em_CA <- c(em_CA,
               future_em_CA) # EM Committee Averaging
    
    # Enfin, on crée/complète le stack de Species Range Change (SRC)
    if(s == scenarios[1])
    {
      SRC_stack <- SRCRaster
    } else
    {
      SRC_stack <- c(SRC_stack,
                     SRCRaster)
    }
  }
  
  names(em_cont) <- c('baseline', scenarios)
  names(em_sd) <- c('baseline', scenarios)
  names(em_bin) <- c('baseline', scenarios)
  names(em_penalised_cont) <- c('baseline', scenarios)
  names(em_penalised_bin) <- c('baseline', scenarios)
  names(em_CA) <- c('baseline', scenarios)
  names(SRC_stack) <- scenarios
  
  cat(paste("\n-- Writing stacks to hdd\n", sep = " "))
  writeRaster(em_cont, paste0("outputs/em_cont_", sp, ".tif"), overwrite = T)
  writeRaster(em_sd, paste0("outputs/em_sd_", sp, ".tif"), overwrite = T)
  writeRaster(em_bin, paste0("outputs/em_bin_", sp, ".tif"), overwrite = T)
  writeRaster(em_penalised_cont, paste0("outputs/em_penalised_cont_", sp, ".tif"), overwrite = T)
  writeRaster(em_penalised_bin, paste0("outputs/em_penalised_bin_", sp, ".tif"), overwrite = T)
  writeRaster(em_CA, paste0("outputs/em_CA_", sp, ".tif"), overwrite = T)
  writeRaster(SRC_stack, paste0("outputs/src_", sp, ".tif"), overwrite = T)

  rm(list = ls()[grep(sp, ls())])
  cat(paste("----", sp, "stack creation finished ----\n", sep = " "))
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
library(terra)

if(!dir.exists("outputs"))
{
  dir.create("outputs", recursive = T)
}
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")


proj_names <- readRDS("data/projection_names.RDS")

stacks <- c("baseline",
            apply(proj_names[, c("gcm_lu", "scenar_clim", "year_lu")], 1, 
                  paste, collapse = "_"))


scenarios <- c("ssp245_2070", "ssp585_2070")

allsp_cutoffs <- readRDS("data/jaccard_cutoffs.RDS")



#### NOTE ####
# This script creates a lot of ensemble projections, but YOU DO NOT NEED TO
# USE ALL OF THEM. Choose only the ones that you need for your specific study.
# Conceptualise your study BEFORE using this script and decide what you want
# to use.


# 7 stacks will be created:  

# - em_cont  # Ensemble model suitabilities
# - em_sd # Ensemble model uncertainty (sd)
# - em_bin  # Ensemble model potential P/A
# - SRC_stack # Species Range Change stack based on em_bin

# - em_penalised_cont # Penalised EM suitabilities
# - em_penalised_bin # Penalised EM potential P/A


# - em_CA # EM Committee Averaging



# Note: a penalised SRC stack is not computed, as it does not make sense:
# Current values are based on a single climate projection
# Future values are based on 3 climate projection
# Hence, the standard deviation in the future is increased by the use of different GCM
# Comparing current values to future values would overestimate the change because of 
# the increased variability between models + GCMs


for (i in 1:nrow(sp_list))
{
  sp <- sp_list$sp[i]
  cat(paste("----", Sys.time(), sp, "stack creation initialised ----\n", sep = " "))
  
  
  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))
  # Récupération des données entrées (notamment pour avoir les coordonnées 
  # des pts de présence)
  input_data <- get_formal_data(model_runs)
  sp_coords <- input_data@coord
  
  
  # Récupération des seuils issus de Jaccard
  tmp <- reshape2::melt(allsp_cutoffs[sp, , , ])
  cur_cutoffs <- tmp$value
  # On va s'assurer que les noms pour les cutoffs correspondent bien aux noms
  # des projections faites par biomod
  names(cur_cutoffs) <- paste(sp, tmp$pa.run, tmp$cv.run, tmp$model, sep = "_")
  
  ##### Partie 1 : baseline #####
  # Cartes de probabilité de présence
  cur_stack <- rast(paste0("models/", sp, "/proj_", stacks[1],
                           "/proj_", stacks[1], "_", sp, ".tif"))
  # Transformation en présence-absence par le seuil de Jaccard
  cur_binary_stack <- bm_BinaryTransformation(cur_stack,
                                              cur_cutoffs[names(cur_stack)])
  
  # Ensemble model: probabilités
  cur_em <- mean(cur_stack) # Ensemble model
  cur_sd <- app(cur_stack, sd) # Uncertainty (standard deviation)
  cur_em_penalised <- cur_em - cur_sd # Penalised probability of presence
  cur_em_penalised[cur_em_penalised < 0] <- 0
  
  # Ensemble model: conversion en présence-absence
  # Ici j'illustre avec l'indice de Jaccard
  # Mais il est facile de reprendre plutôt le TSS ou ROC 
  # en allant chercher les évaluations dans l'objet d'EM de biomod
  
  
  obs_data <- input_data@data.species
  obs_data[is.na(obs_data)] <- 0 # Pour les présences pseudoabs,
  # On remplace les pseudoabs par des absences
  # Notez bien que c'est une hypothèse très forte (très fausse), mais 
  # c'est la seule possibilité en presence pseudoabs
  # A éviter si vous pouvez vous en passer. Soyez inventifs sur votre 
  # utilisation des sorties continues des modèles !
  # Cf. e.g. la méthode que j'ai utilisé dans le rapport SDMs Corse
  # https://borisleroy.com/sdms-pna-corse/
  
  # On extrait les probas aux coordonnées de présence
  # Notez qu'on le fait sur tous les points de présence (contrairement au script
  # Jaccard où on ne le faisait que sur les données d'évaluation)
  # Ce n'est pas un problème ici car **on n'évalue pas le modèle** : 
  # on recherche simplement le seuil optimal
  pred_data <- extract(cur_em, sp_coords)
  jaccard_test <- NULL
  for(cutoff in seq(0, 1000, by = 1))
  {
    pred.pa <- pred_data$mean
    pred.pa[pred.pa < cutoff] <- 0
    pred.pa[pred.pa >= cutoff] <- 1
    TP <- length(which(obs_data == 1 & pred.pa == 1))
    FN <- length(which(obs_data == 1 & pred.pa == 0))
    FP <- length(which(obs_data == 0 & pred.pa == 1))
    jaccard <- TP / (TP + FP + FN)
    jaccard_test <- rbind.data.frame(jaccard_test,
                                     data.frame(cutoff = cutoff,
                                                TP = TP,
                                                FN = FN,
                                                FP = FP,
                                                jaccard = jaccard))
  }
  em_cutoff <- mean(jaccard_test$cutoff[
    which(jaccard_test$jaccard == max(jaccard_test$jaccard))
  ])
  
  # Ensemble model transformé en présence-absence
  cur_em_binary <- bm_BinaryTransformation(cur_em, em_cutoff)
  # Transformation en présence-absence des probas pénalisées par l'incertitude
  cur_em_binary_penalised <- bm_BinaryTransformation(cur_em_penalised, 
                                                     em_cutoff)
  # Committee averaging
  cur_em_CA <- sum(cur_binary_stack)
  cur_em_CA[cur_em_CA < (nlyr(cur_binary_stack) / 2)] <- 0
  cur_em_CA[cur_em_CA >= (nlyr(cur_binary_stack) / 2)] <- 1
  
  
  # On prépare nos stacks qui contiendront tous les rasters finaux ici
  em_cont <- cur_em # Ensemble model suitabilities
  em_sd <- cur_sd # Ensemble model uncertainty (sd)
  em_bin <- cur_em_binary # Ensemble model potential P/A
  em_penalised_cont <- cur_em_penalised # Penalised EM suitabilities
  em_penalised_bin <- cur_em_binary_penalised # Penalised EM potential P/A
  em_CA <- cur_em_CA # EM Committee Averaging
  
  
  ##### Partie 2 : scénarios futurs #####
  # La difficulté ici est qu'il faut réunir les projections des différents GCMs,
  # stockés par biomod dans différents objets 
  for(s in scenarios)
  {
    cat(paste("\n-- Calculations for scenario:", s, "\n", sep = " "))
    # Il faut tout d'abord réunir les résultats de tous les
    # GCMs de notre scénario dans un seul et unique stack
    subStacks <- stacks[grep(s, stacks)]
    for (s1 in subStacks)
    {
      if(s1 == subStacks[1])
      {
        future_stack <- rast(paste("models/", sp, "/proj_", s1, "/", 
                                   "proj_", s1, "_", sp, ".tif", sep =""))
        order_list <- names(future_stack) # On récupère les noms que biomod a
        # donné, afin de garder le lien entre modèles et cutoffs à appliquer
        
        # Ensuite on renomme les couches en ajoutant le nom du GCM au début
        names(future_stack) <-  paste(s1, names(future_stack), sep = "_")
      } else
      {
        tmp <- rast(paste("models/", sp, "/proj_", s1, "/",
                          "proj_", s1, "_", sp, ".tif", sep =""))
        names(tmp) <- paste(s1, names(tmp), sep = "_")
        future_stack <- c(future_stack,
                          tmp)
      }
    }
    ##### Attention: si votre protocole est trop lourd (trop de modèles, 
    # trop de répétitions, trop de GCMs, taille de la zone d'étude trop grande 
    # ou résolution trop fine), les étapes suivantes pourraient être
    # trop lourdes pour un pc moyen. Dans ce cas il faudra
    # se "contenter" de faire la moyenne entre les modèles 
    # d'ensemble projetés par biomod sur les différents GCM
    
    # Cartes de présence-absence pour chaque modèle du futur
    future_binary_stack <- bm_BinaryTransformation(future_stack, 
                                                   rep(cur_cutoffs[order_list], 
                                                       length(subStacks)))
    
    # Ensemble model: probabilités
    future_em <- mean(future_stack) # Ensemble model
    future_sd <- app(future_stack, sd) # Uncertainty (standard deviation)
    future_em_penalised <- future_em - future_sd # Penalised probability of presence
    future_em_penalised[future_em_penalised < 0] <- 0
    
    # Ensemble model transformé en présence-absence
    future_em_binary <- bm_BinaryTransformation(future_em, em_cutoff)
    
    # Transformation en présence-absence des probas pénalisées par l'incertitude
    future_em_binary_penalised <- bm_BinaryTransformation(future_em_penalised, 
                                                          em_cutoff)
    
    # Committee averaging
    future_em_CA <- sum(future_binary_stack)
    future_em_CA[future_em_CA < (nlyr(future_binary_stack) / 2)] <- 0
    future_em_CA[future_em_CA >= (nlyr(future_binary_stack) / 2)] <- 1
    
    # Range size change (basé sur les rasters NON PENALISES)
    # WARNING: ne pas utiliser les rasters pénalisés ici, car
    # l'incertitude est plus grande dans le futur à cause de 
    # l'utilisation de différents GCMs
    # De toute façon à ce stade j'espère vous avoir convaincu de favoriser 
    # la proba plutôt que les stacks binaires !!!
    SRCRaster <- future_em_binary - 2 * cur_em_binary
    # -2 : aire perdue
    # -1 : aire stable
    # 0 : aire inoccupée dans les 2 périodes
    # 1 : aire nouvelle
    
    # tmp <- SRCRaster
    # tmp[tmp != -2] <- 0
    # tmp[tmp == -2] <- 1
    # global(tmp * cellSize(tmp), fun = "sum", na.rm = TRUE)
    

    # On complète nos stacks avec les résults de l'EM pour le scénario actuel
    em_cont <- c(em_cont,
                 future_em) # Ensemble model suitabilities
    em_sd <- c(em_sd,
               future_sd) # Ensemble model uncertainty (sd)
    em_bin <- c(em_bin,
                future_em_binary) # Ensemble model potential P/A
    em_penalised_cont <- c(em_penalised_cont,
                           future_em_penalised) # Penalised EM suitabilities
    em_penalised_bin <- c(em_penalised_bin,
                          future_em_binary_penalised) # Penalised EM potential P/A
    em_CA <- c(em_CA,
               future_em_CA) # EM Committee Averaging
    
    # Enfin, on crée/complète le stack de Species Range Change (SRC)
    if(s == scenarios[1])
    {
      SRC_stack <- SRCRaster
    } else
    {
      SRC_stack <- c(SRC_stack,
                     SRCRaster)
    }
  }
  
  names(em_cont) <- c('baseline', scenarios)
  names(em_sd) <- c('baseline', scenarios)
  names(em_bin) <- c('baseline', scenarios)
  names(em_penalised_cont) <- c('baseline', scenarios)
  names(em_penalised_bin) <- c('baseline', scenarios)
  names(em_CA) <- c('baseline', scenarios)
  names(SRC_stack) <- scenarios
  
  cat(paste("\n-- Writing stacks to hdd\n", sep = " "))
  writeRaster(em_cont, paste0("outputs/em_cont_", sp, ".tif"), overwrite = T)
  writeRaster(em_sd, paste0("outputs/em_sd_", sp, ".tif"), overwrite = T)
  writeRaster(em_bin, paste0("outputs/em_bin_", sp, ".tif"), overwrite = T)
  writeRaster(em_penalised_cont, paste0("outputs/em_penalised_cont_", sp, ".tif"), overwrite = T)
  writeRaster(em_penalised_bin, paste0("outputs/em_penalised_bin_", sp, ".tif"), overwrite = T)
  writeRaster(em_CA, paste0("outputs/em_CA_", sp, ".tif"), overwrite = T)
  writeRaster(SRC_stack, paste0("outputs/src_", sp, ".tif"), overwrite = T)

  rm(list = ls()[grep(sp, ls())])
  cat(paste("----", sp, "stack creation finished ----\n", sep = " "))
>>>>>>> 5768564 (2nd commit)
}