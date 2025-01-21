
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
# Using multi-dimension arrays for md eval :
# Array 2D : md vs sp for a given run of background
# Array 3D : md vs sp for all runs of background for a given CV run
# Array 4D : md vs sp for all runs of background and for all runs of CV

# To extract info from multiD array:
# d["sp1",,,"cv1 "] --> for species 1, for all background runs, for 1rst CV run.

# Load libraries
library(biomod2)

# Load species list
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")


# On charge arbitrairement la première espèce avec PA pour récupérer
# automatiquement les paramètres de notre modélisation : 
# algorithmes, runs de CV, runs de PA.
# Nom de l'espèce
sp <- sp_list$sp[which(sp_list$pa.generation == "yes")][1]

# Objet avec les modèles calibrés
model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))

# On crée un array dans lequel seront stockés les évaluations Jaccard de tous les modèles 
# Toutes les spp doivent avoir les même modèles et les mêmes runs, sinon il faut rentrer qq infos à la main (voir les commentaires et l'exemple commenté ci-dessous)
allsp_jaccard <- array(
  
  # L'array va avoir 4 dimensions :
  dim = c(
    # Les espèces
    nrow(sp_list),
    # Les algorithmes de modélisation
    length(unique(model_runs@models.evaluation@val$algo)), # Si on a des modèles différents pour chaque espèce, on rentre le nombre de modèle total à la main ici 
    # Les runs de cross-validation      
    length(unique(model_runs@models.evaluation@val$run)),
    # Les runs de pseudo-absences
    length(unique(model_runs@models.evaluation@val$PA)) + 1),
  
  # On donne les noms pour chaque de ces dimensions:
  dimnames = list(
    # Espèces
    species = sp_list$sp,
    # Algorithmes
    model = unique(model_runs@models.evaluation@val$algo), # Pareil, ici si on a différents modèles selon les spp, on les met ici à la main
    # Runs de CV              
    cv.run = unique(model_runs@models.evaluation@val$run),
    # Runs de PA
    pa.run = c("allData", # Quand on est en présence-absence ça s'appelle 
               # toujours 'allData' dans biomod2
               unique(model_runs@models.evaluation@val$PA))))

# Un autre facçon de faire la même chose, en plus simple, mais qui demande plus 
# de remplissage à la main, serait de repréciser ici nos paramètres de
# modélisation:
# models <-  c("GLM", "GAM", "ANN", "GBM", "FDA", "RF", "MAXNET", "XGBOOST")
# runs_CV <- 2
# runs_PA <- 2
# allsp_jaccard <- array(
#   # L'array va avoir 4 dimensions :
#   dim = c(
#     # Les espèces
#     nrow(sp_list),
#     # Les algorithmes de modélisation
#     length(models),
#     # Les runs de cross-validation
#     runs_CV,
#     # Les runs de pseudo-absences
#     runs_PA + 1),
# 
#   # On donne les noms pour chaque de ces dimensions:
#   dimnames = list(
#     # Espèces
#     species = sp_list$sp,
#     # Algorithmes
#     model = models,
#     # Runs de CV
#     cv.run = paste0("RUN", 1:runs_CV),
#     # Runs de PA
#     pa.run = c("allData", # Quand on est en présence-absence ça s'appelle
#                # toujours 'allData' dans biomod2
#                paste0("PA", 1:runs_PA))))



# On crée également un array dans lequel seront stockés les seuils optimaux de 
# conversion en présence-absence, afin de pouvoir les utiliser dans le 
# modèle d'ensemble
allsp_cutoffs <- allsp_jaccard


# Et on prépare une liste qui va contenir le détail des calculs au cas où on veut y revenir plus tard
jaccard_list <- list()

# Note: la boucle n'est pas optimisée pour avoir un code rapide et le plus court
# possible, mais ce n'est pas un problème ici car le calcul est généralement
# extrêmement rapide.
# L'objectif ici et de décomposer le code pour que vous compreniez bien le
# processus de calcul des indices.

for (i in 1:nrow(sp_list)) {
  
  sp <- sp_list$sp[i]
  cat(paste("----", Sys.time(), sp, " evaluation initialised ----\n", sep = " "))
  
  # Récupération de tout ce qu'il faut pour évaluer nos modèles : 
  # occurrences, prédictions des modèles, et répartition des occurrence en
  # calibration et évaluation

  # Récupération des objets avec les modèles calibrés
  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))
  
  # Récupération des données d'occurrence 
  input_data <- get_formal_data(model_runs)
  
  # Récupération des lignes de calibration et évaluation
  calib_lines <- get_calib_lines(model_runs) # Determine pour chaque point et pour chauque calib run s'il a été utilisé (TRUE) ou non (FALSE)
  
  # Récupération des prédictions des modèles (/!\ si pred = 2, proba = 0.002, pred = 1000, proba = 1)
  model_preds <- get_predictions(model_runs)

  # On fait maintenant une boucle pour travailler sur chaque run de PA
  for (pa in unique(model_preds$PA))
  {
    
    # Ensuite, une boucle pour travailler sur chaque run de CV
    for (cv in unique(model_preds$run)) {
      
      # Etape 1 : préparer les observations qui servent à l'évaluation
      # On va isoler les données d'occurrence qui doivent servir à l'évaluation
      # = données d'occurrences  qui n'ont pas servi à calibrer le modèle pour 
      # le run actuel de PA et de CV
      obs_data <- input_data@data.species[
        which(!calib_lines[, paste0("_", pa, "_", cv)])
      ]
      
      # /!\ Transformation pseudo-absences en absences (pas bien !)
      # On ne peut PAS calculer l'indice de Jaccard sur les présences seules 
      # Ici on a des présences et des abscences (0)
      # Mais si on voulait tester cet indice avec un jeu de donnée en presence seule (ie presence + background = NA) (pas bien !), il faudrait remplacer les background (NA) en abscence (0):
      obs_data[is.na(obs_data)] <- 0 
      
      # Etape 2 : préparer les prédictions des modèles à évaluer
      # On récupère ici les lignes qui ont servi à la calibration des modèles
      cur_calib_lines <- calib_lines[
        , paste0("_", pa, "_", cv)
      ]
      
      # Il peut y avoir des NAs dans ce vecteur, que l'on élimine ici,
      # ce qui est nécessaire car biomod aura éliminé en interne ces NAs
      # également (ce qui peut alors causer un mismatch dans notre analyse)
      cur_calib_lines <- cur_calib_lines[which(!is.na(cur_calib_lines))]
      
      # Eventuellement print s'il y a des NA (using any()):
     # print(any(is.na(cur_calib_lines)))
      
      # Une fois que l'on sait quelles sont les lignes de calibration, on peut
      # automatiquement retrouver les numéros (indices) des lignes d'évaluation
      cur_eval <- which(!cur_calib_lines)
      
      for (md in unique(model_preds$algo))
      {
        
        # On évalue uniquement les lignes qui n'ont pas servi à la calibration :
        cur_preds <- model_preds[
          which(model_preds$points %in% cur_eval & # Lignes d'évaluation uniquement 
                  model_preds$algo == md & # Algorithme actuel
                  model_preds$run == cv & # run cv actuel
                  model_preds$PA == pa), ] # run PA actuel
        
        # S'il y a des NAs ou pas de lignes dans cur_preds,
        # le modèle a échoué donc on ne calcule pas jaccard
        if(!any(is.na(cur_preds)) & !(nrow(cur_preds) == 0)) 
        {
          
          # Calcul de l'indice de Jaccard pour tous les seuils entre 0 et 1
          jaccard_test <- NULL
          
          for(cutoff in seq(0, 1000, by = 1))
          {
            
            # On convertit les probas en binaire 1/0 en fonction du cutoff
            pred_pa <- cur_preds$pred
            pred_pa[pred_pa < cutoff] <- 0
            pred_pa[pred_pa >= cutoff] <- 1
            
            # On calcule les éléments de la matrice de confusion
            # True Positives
            TP <- length(which(obs_data == 1 & pred_pa == 1)) 
            
            # False Negatives
            FN <- length(which(obs_data == 1 & pred_pa == 0))
            
            # False Positives
            FP <- length(which(obs_data == 0 & pred_pa == 1))
            
            # TN <- length(which(obs_data == 0 & pred_pa == 0))
            jaccard <- TP / (TP + FP + FN)
            jaccard_test <- rbind.data.frame(jaccard_test,
                                             data.frame(cutoff = cutoff,
                                                        TP = TP,
                                                        FN = FN,
                                                        FP = FP,
                                                        jaccard = jaccard))
          }
          jaccard_list[[sp]][[paste0(pa, "_", cv, "_", md)]] <- jaccard_test
          
          # Le seuil retenu est celui qui a le plus fort indice de Jaccard
          # On prend la moyenne du meilleur seuil, si celui-ci est atteint plusieurs fois
          allsp_cutoffs[sp, md, cv, pa] <- mean(jaccard_test$cutoff[
            which(jaccard_test$jaccard == max(jaccard_test$jaccard))])
          
          # On extrait la valeur de jaccard au meilleur seuil moyen
          allsp_jaccard[sp, md, cv, pa] <- jaccard_test$jaccard[
            which(jaccard_test$cutoff == round(allsp_cutoffs[sp, md, cv, pa]))]
        } else
        {
          jaccard_list[[sp]][[paste0(pa, "_", cv, "_", md)]] <- NA
          
          allsp_cutoffs[sp, md, cv, pa] <- NA 
          
          allsp_jaccard[sp, md, cv, pa] <- NA
        }
      }
    }
  }
}

# On sauvegarde les cutoffs pour les réutiliser dans le modèle d'ensemble si besoin
saveRDS(allsp_cutoffs, file = "./data/jaccard_cutoffs.RDS")

# Et on sauvegarde les évaluations de Jaccard
saveRDS(allsp_jaccard, file = "./data/jaccard_evals.RDS")
saveRDS(jaccard_list, file = "./data/jaccard_tests.RDS")

##### Partie 2 : Graphiques #####
allsp_jaccard <- readRDS("./data/jaccard_evals.RDS")
library(ggplot2)

ggjaccard <- reshape2::melt(allsp_jaccard)

ggplot(ggjaccard, aes(x = model, y = value, col = pa.run, shape = cv.run)) +
  geom_point() + facet_wrap (~ species)



