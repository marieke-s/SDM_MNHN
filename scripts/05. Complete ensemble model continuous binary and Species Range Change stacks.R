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
library(terra)

# Create outputs directory if it does not exist
if(!dir.exists("outputs"))
{
  dir.create("outputs", recursive = T)
}

# Retrieve species list
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")

# Retrieve projection names
proj_names <- readRDS("data/projection_names.RDS")

# Create stack names
stacks <- c("baseline",
            apply(proj_names[, c("gcm_lu", "scenar_clim", "year_lu")], 1, 
                  paste, collapse = "_"))

# On a 2 scénarios, la question qui se pose est : comment on va combiner les projections ? 
# Ici on met dans le nom le scenario et la période temporelle (car il peut en y avoir plusieurs)
scenarios <- c("ssp245_2070", "ssp585_2070")

# Load Jaccard cutoffs pour la binarisation des projections
allsp_cutoffs <- readRDS("data/jaccard_cutoffs.RDS")



  ##### NOTE #####
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
  
  # Load models
  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))
  
  # Récupération des données entrées (notamment pour avoir les coordonnées 
  # des pts de présence)
  input_data <- get_formal_data(model_runs)
  sp_coords <- input_data@coord
  
  # Récupération des seuils issus de Jaccard (meilleur seuil, cf. script 02b.)
  tmp <- reshape2::melt(allsp_cutoffs[sp, , , ]) # obj temporaire
  cur_cutoffs <- tmp$value
  
  # On va s'assurer que les noms pour les cutoffs correspondent bien aux noms
  # des projections faites par biomod
  names(cur_cutoffs) <- paste(sp, tmp$pa.run, tmp$cv.run, tmp$model, sep = "_")
  
  
  ##### Partie 1 : baseline #####
  
  # Cartes de probabilité de présence
  cur_stack <- rast(paste0("models/", sp, "/proj_", stacks[1],
                           "/proj_", stacks[1], "_", sp, ".tif")) # on récupère les probabilités individuelles de tous les modèles pour le baseline
  
  # ! Attention ! 
  # Certains modèles ont des préditcions qui sorte du range théorique 0-1 : eg. le modèle MARS produit des valeurs négatives !
  # Attention : si on observe des valeurs négatives il faut regarder les courbes de réponses et se poser la question de la validité des résultats. Ce bornage peut être dû a une extrapolation trop forte mais ça peut aussi être un problème de du modèle --> choisir si on élimine le modèle MARS (doute sur l'implémentation du modèle par biomod). 
  # Attention : la correction qu'on applique ci-dessous arrive peut être déjà trop tard dans le pipeline : réfléchir à quand elle devrait être appliquée.
  
  # Pour checker ces valeurs : 
  # a <- values(cur_stack)
  # apply(a,2, min, na.rm = TRUE)
  # apply(a,2, max, na.rm = TRUE)
  # On corrige donc les valeurs pour les mettre à 0

  # Bornage des proba entre 0 et 1000
  # Ici on choisit de garder le modèle sans analyse approfondie du problème. 
  # Transformation des valeurs négatives en 0
  cur_stack[cur_stack < 0] <- 0 
  cur_stack[cur_stack > 1000] <- 1000
  
  # Transformation en présence-absence par le seuil de Jaccard (binarisation) avec bm_BinaryTransformation()
  cur_binary_stack <- bm_BinaryTransformation(data = cur_stack,
                                              threshold = cur_cutoffs[names(cur_stack)]) # on ordonne les valeurs de cutoffs selon l'ordre des modèles dans le stack
  
  # Optionel : plot
  # plot(cur_binary_stack)
  
  # Ensemble model: probabilités
  cur_em <- mean(cur_stack) # Ensemble model (moyenne des probas des modèles)
  cur_sd <- app(cur_stack, sd) # Uncertainty (standard deviation)
  cur_em_penalised <- cur_em - cur_sd # Penalised probability of presence (proba - sd)
  cur_em_penalised[cur_em_penalised < 0] <- 0 # Bornage à 0 car on ne veut pas que de proba négative (même si sd élevé)
  
  # Optionel : plot both cur_em_penalised and cur_em side to side
  # plot(c(plot(cur_em), plot(cur_em_penalised)))
  # On observe par ex qu'à l'est, les habitats fav sont pénalisés --> carte restrictive produite --> à utiliser dans le cadre de prise de décision. 

  # Ensemble model: conversion en présence-absence
  # Ici j'illustre avec l'indice de Jaccard
  # Mais il est facile de reprendre plutôt le TSS ou ROC 
  # en allant chercher les évaluations dans l'objet d'EM de biomod
  
  # Récupération du vecteur de P-A
  obs_data <- input_data@data.species
  
  # Si on a des présences et des pseudoabs,
  # On remplace les pseudoabs par des absences
  # Notez bien que c'est une hypothèse très forte (très fausse), mais 
  # c'est la seule possibilité en presence pseudoabs
  # A éviter si vous pouvez vous en passer. Soyez inventifs sur votre 
  # utilisation des sorties continues des modèles !
  # Cf. e.g. la méthode que j'ai utilisé dans le rapport SDMs Corse
  # https://borisleroy.com/sdms-pna-corse/
  obs_data[is.na(obs_data)] <- 0 
  
  # On extrait les probas aux coordonnées de présence
  # Notez qu'on le fait sur tous les points de présence (contrairement au script
  # Jaccard où on ne le faisait que sur les données d'évaluation)
  # Ce n'est pas un problème ici car **on n'évalue pas le modèle** : 
  # on recherche simplement le seuil optimal
  pred_data <- extract(cur_em, sp_coords) # valeur de probabilités moyennes 
  jaccard_test <- NULL
  
  # Ici on calcul le jaccard sur le md d'esemble sur le jeu de donnée qu'il a déjà vu mais c'est pas pour en faire l'évaluation, c'est uniquement pour trouver le meilleur seuil de binarisation.
  
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
  
  # On récupère les valeurs de l'indice de jaccard en fct du cutoff
  em_cutoff <- mean(jaccard_test$cutoff[
    which(jaccard_test$jaccard == max(jaccard_test$jaccard)) # on récupère le cutoff pour lequel on a le meilleur jaccard.
  ])
  
  # Ensemble model transformé en présence-absence
  cur_em_binary <- bm_BinaryTransformation(cur_em, em_cutoff)
  
  # Transformation en présence-absence des probas pénalisées par l'incertitude
  cur_em_binary_penalised <- bm_BinaryTransformation(cur_em_penalised, 
                                                     em_cutoff)
  # Committee averaging
  cur_em_CA <- sum(cur_binary_stack)
  cur_em_CA[cur_em_CA < (nlyr(cur_binary_stack) / 2)] <- 0 # presence à partir de 50% des modèles ont une probabilité au dessus de cutoff, ie = 1)
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
    subStacks <- stacks[grep(s, stacks)] # substack c'est les noms des stacks des modèles filtrés pour nos scénarios
    
    for (s1 in subStacks)
    {
      if(s1 == subStacks[1])
      {
        # Lire la projection future qui correspond
        future_stack <- rast(paste("models/", sp, "/proj_", s1, "/", 
                                   "proj_", s1, "_", sp, ".tif", sep =""))
        
        # On récupère les noms que biomod a
        # donné, afin de garder le lien entre modèles et cutoffs à appliquer
        # pour faire la correspondance avec les modèles en binaire
        order_list <- names(future_stack) 
        
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
    # --> Avec terra c'est peut-être plus le cas.
    
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
    
    # Species range change (difference present - futur) :
    # Range size change (basé sur les rasters NON PENALISES)
    # WARNING: ne pas utiliser les rasters pénalisés ici, car
    # l'incertitude est plus grande dans le futur à cause de 
    # l'utilisation de différents GCMs
    # De toute façon à ce stade j'espère vous avoir convaincu de favoriser 
    # la proba plutôt que les stacks binaires !!!
    SRCRaster <- future_em_binary - 2 * cur_em_binary
    
    # Interprétation:
    # -2 : aire perdue
    # -1 : aire stable
    # 0 : aire inoccupée dans les 2 périodes
    # 1 : aire nouvelle
    
    # Binarisation 
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
    if(s == scenarios[1]) # création du stack pour la 1ere
    {
      SRC_stack <- SRCRaster
    } else                # complétion du stack pour les suivantes
    {
      SRC_stack <- c(SRC_stack,
                     SRCRaster)
    }
  }
  
  # Renommage de toutes les couches
  names(em_cont) <- c('baseline', scenarios)
  names(em_sd) <- c('baseline', scenarios)
  names(em_bin) <- c('baseline', scenarios)
  names(em_penalised_cont) <- c('baseline', scenarios)
  names(em_penalised_bin) <- c('baseline', scenarios)
  names(em_CA) <- c('baseline', scenarios)
  names(SRC_stack) <- scenarios
  
  # Sauvegarde
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
}
