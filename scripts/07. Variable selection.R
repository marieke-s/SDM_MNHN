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

# Préambule :
# Avant d'utiliser ce script il faut au maximum sélectionner vos variables
# sur la base d'un raisonnement biologique. Ainsi, il ne faudrait 
# théoriquement jamais le lancer sur l'ensemble des variables comme on va le
# faire dans l'exemple ici. Il faudrait plutôt avoir sélectionné a priori
# quelles variables climatiques, quelles variables d'habitat, etc. sont 
# importantes pour vos espèces ou votre groupe d'espèces, en étayant votre
# raisonnement avec de la bibliographie.

# Load libraries
library(tidyr)  
library(terra)
library(biomod2)
library(virtualspecies)
library(Rarity) # Il existe des packages plus à jour
library(ggplot2)
library(usdm)

# Load baseline data
baseline <- rast("data/baseline.tif")

# Load species list
sp_list <- read.csv("./data_cours/species_list.csv", sep = ";")

# A noter : la sélection de variable doit avant tout se faire sur la base
# des connaissances écophysiologiques de l'espèce - une sélection automatique
# comme dans ce script ne doit intervenir que dans un second temps sur les 
# variables pré-sélectionnées. Pour voir un exemple détaillé de sélection 
# sur la base d'hypothèses sur les espèces, cf. le rapport SDMs Corse 
# https://borisleroy.com/sdms-pna-corse/

# Entrer la liste des modèles sur lesquels travailler
# (idéalement la même que les modèles que vous visez au final)
models <- c('GLM','GAM', 'ANN', 'RF', 'FDA', 'MAXNET', 'XGBOOST')

# Recherche de groupes de variables intercorrélées
png("graphiques/collinearity_groups.png") # sépare les variables en groupes de variables intercorrelées 
groups <- removeCollinearity(baseline, plot = T,
                             multicollinearity.cutoff = 0.7, # Dormann et al. 2007 : justifie le seuil de 0.7 # Attention, si vous avez un raisonnement fort pour garder 2 variables faites le malgrès le seuil.
                             method = "spearman") 

dev.off()

# Evaluation des variance inflation factor
# Normalement le jeu final de variables ne devrait pas dépaser 10 en VIF
# Notez que les fonctions comme vifstep agissent de manière brutale et peuvent
# éliminer les bonnes variables du dataset. Mieux vaut prioriser le savoir
# biologique et choisir les variables qui ont le plus de sens / sur lesquelles
# vous pouvez poser une hypothèse pour vos espèces. 
vifs <- usdm::vif(baseline) # Toutes les variables

# Retenir une seule par groupe de variables intercorrélées
vifs <- usdm::vif(baseline[[c("bio3", "bio2", "bio4", 
                              "bio5", "bio12", "forests")]])

# Echanger les variables entre elles dans chaque groupe permet de réduire le vif
vifs <- usdm::vif(baseline[[c("bio3", "bio2", "bio19", 
                              "bio5", "bio12", "forests")]])

# Sélection automatique de variables
vifs <- usdm::vifstep(baseline[[c("bio6", "bio2", "bio19", 
                                  "bio5", "bio12", "forests")]]) # on fait un vifstep sur un jeu réduit de variable que l'on présuppose importante pour l'espèce.

# Peut-être faut-il plutôt explorer le VIF sur les données de présence-absence
# plutôt que sur tout le raster ?

# On met dans un objet tous les groupes qui ont plusieurs variables
sel_groups <- groups[which(sapply(groups, length) > 1)]

# Fin VIF #

# On prépare un tableau qui va contenir les variables sélectionnées pour chaque 
# groupe
kept_vars <- data.frame(species = sp_list$sp)

##### Etape 1 : sélection d'une variable par groupe ####

# Ici je fais une sélection arbitraire manuelle
# Choisissez toujours la variable qui a le plus de sens biologique pour votre
# espèce.
kept_vars$group1 <- "bio6"
kept_vars$group2 <- "bio7"
kept_vars$group3 <- "bio5"
kept_vars$group4 <- "bio16"
kept_vars$group5 <- "bio17"

# Attention, ici on choisit les mêmes variables pour toutes les espèces, c'est possible de mettre différentes variables pour chaque espèce.

# On récupère maintenant toutes les variables pas corrélées aux autres
lone_vars <- unlist(groups[which(sapply(groups, length) == 1)])

# Et on va les rajouter aux variables sélectionnées dans chaque groupe de 
# variables intercorrélées
if(nrow(kept_vars) == 1)
{
  kept_vars[, lone_vars] <- lone_vars
  vars_per_sp <- kept_vars
} else
{
  vars_per_sp <- data.frame(kept_vars,
                            sapply(lone_vars, rep,  nrow(kept_vars)))
}

# Save the variables
saveRDS(vars_per_sp, file = paste0("data/vars_per_sp.RDS"))


##### Etape 2: calibration d'un grand modèle avec toutes les vars sélectionnées #####
vars_per_sp <- readRDS(paste0("data/vars_per_sp.RDS"))

for (i in 1:length(sp_list$sp))
{
  sp <- sp_list$sp[i]
  P_points <- readRDS(paste0("./data/occurrences_", sp))  
  sp_env_stack <- baseline[[
    vars_per_sp[vars_per_sp$species == sp, 
                             2:ncol(vars_per_sp)]]]
  
  # 1. Biomod2 doit-il générer des pseudoabsences ?
  if(sp_list$pa.generation[i] == "yes")
  {
    nb_PA <- 1000 # Générer autant de pseudoabsences
    # que de présences
    runs_PA <- 2 # Nombre de répétitions de pseudo-absences
    runs_CV <- 2 # Nombre de runs de CV
    PA_strategy <- "random" # Stratégie de génération des pseudo-absences
  } else
  {
    runs_PA <- 0
    nb_PA <- 0
    runs_CV <- 4 # On monte le nombre de runs de CV dans le cas des # Si on est en P-A il faut augmenter le nb de CV pour augmenter les tests --> on teste les variables donc il faut plus de run.
    # P/As pour avoir plus de répétitions
    PA_strategy <- NULL
  }
  
  
  # 2. Formatage des occurrences pour biomod
  coorxy <- P_points[, c("x", "y")]
  P_points <- P_points[, "Observed"]

  # Creation du dossier de stockage des résultats
  if(!dir.exists(paste0("var_selection/", sp)))
  {
    dir.create(paste0("var_selection/", sp), recursive = T)
  }

  # Formatage des données pour biomod
  run_data <- BIOMOD_FormatingData(resp.name = sp,
                                   resp.var = P_points, 
                                   expl.var = sp_env_stack[[
                                     sample(names(sp_env_stack), 
                                            nlyr(sp_env_stack))]], 
                                   dir.name = "var_selection",
                                   resp.xy = coorxy, 
                                   PA.nb.rep = runs_PA,
                                   PA.nb.absences = nb_PA, 
                                   PA.strategy = PA_strategy)
  
  # Sauvegarde des données formatées pour biomod
  saveRDS(run_data, file = paste0("var_selection/", 
                                  sp, "/run_data.RDS"))
  
  ##### Important : ajoutez ici la paramétrisation de vos modèles #####
  # Il faut mettre les mêmes paramètres que ce qui va être utilisé par la suite
  # model_parameters <- bm_ModelingOptions()
  
  # Calibration des modèles (ici paramètres un peu par défault)
  model_runs <- BIOMOD_Modeling(
    bm.format = run_data, 
    modeling.id = "1", 
    models =  models,
    CV.strategy = "random", # Stratégie de validation croisée
    CV.nb.rep = runs_CV, # Nombre de runs de validation croisée
    CV.perc = 0.8, # ratio gardé pour la calibration dans la CV
    CV.do.full.models = FALSE, # Faut-il faire des modèles finaux avec 100% des
    # données ?
    OPT.data.type = "binary", # Type de données en input
    OPT.strategy = "bigboss", # Quels paramètres pour les modèles ?
    OPT.user.val = NULL, # Vos paramètres customisés pour les modèles ici
    OPT.user.base = "default", # Quels paramètres par défaut pour les modèles
    # pour les paramètres que vous n'avez pas ajustés ?
    # bm.options = , # Options de modélisation à mettre dans cet argument
    weights = NULL, 
    prevalence = 0.5, 
    var.import = 10, 
    nb.cpu = 4, 
    do.progress = TRUE)

  saveRDS(model_runs, file = paste0("var_selection/", sp, 
                                    "/model_runs.RDS"))
  
}


##### Etape 3 : sélection des variables à importance significative #####
sel_vars <- list()

for (i in 1:length(sp_list$sp))
{
  sp <- sp_list$sp[i]
  
  model_runs <- readRDS(paste0("var_selection/", sp, 
                               "/model_runs.RDS"))
  
  gg_varimp <- get_variables_importance(model_runs)
  
  colnames(gg_varimp) <- c("id", 
                           "PA.Run",
                           "CV.Run", "Model", 
                           "Variable", 
                           "VI.run",
                           "Variable.importance")
  
  # Visualisation de l'importance des variables
  
  # Ordonne les variable selon l'importance mediane
  gg_varimp$Variable <- reorder(gg_varimp$Variable,
                                gg_varimp$Variable.importance,
                                median,
                                na.rm=TRUE)
  
  # Plot variable importance for all variables and all models
  p <-  ggplot(gg_varimp, aes(y = Variable, x = Variable.importance)) +
    geom_boxplot() + 
    geom_jitter(alpha = .5, aes(col = Model)) +
    theme_bw() + 
    ggtitle(sp) + 
    scale_color_brewer(palette = "Set2")
  
  # Plot per model
  p2 <-  ggplot(gg_varimp, aes(y = Variable, x = Variable.importance)) +
    geom_boxplot() + 
    geom_jitter(alpha = .5, aes(col = Model)) +
    facet_wrap(~ Model) +
    theme_bw() + 
    ggtitle(sp) + 
    scale_color_brewer(palette = "Set2")
 
  # Sauvegarde des plots
  png(paste0("graphiques/variable_importance_", sp, ".png"))
  print(p)
  dev.off
  
  # Selection des variables
  
  # Ordrer variables par ordre d'importance
  gg_varimp$Variable <- factor(gg_varimp$Variable, 
                               levels = rev(levels(gg_varimp$Variable)))

  # On ne garde que les variables importantes à plus de x% pour au moins 50% des modèles
  # Alternative : garder les X meilleures variables
  quantile.50 <- aggregate(Variable.importance ~ Variable, 
                           data = gg_varimp, FUN = median)
  
  # Evaluation des interactions entre variables (pour éviter les interactions négatives entre les variables)
  
  # Convert table from long to wide format
  var_interactions <- spread(gg_varimp, Variable, Variable.importance)
  
  # Plot correlation matrix
  png(paste0("graphiques/variable_interactions_", sp, ".png"),
      width = 900, height = 900)
  
  corPlot(var_interactions[, -c(1:5)], method = "pearson") # col 1-5 : noms de md etc.
  
  dev.off()  
  
  if(any(cor(var_interactions[, -c(1:5)], use = "na.or.complete") <= -0.4))
  {
    warning(paste0("Negative interactions among variables have been detected for
                   species ", sp, ".\nYou should check the correlation plots among
                   variable importances."))
    # Retirer les variables de manière itérative: retirer la variable avec 
    # interaction négative à l'importance la plus faible, puis recalculer les 
    # corrélations, et ainsi de suite.
  }
    
  sel_vars[[sp]] <- as.character(quantile.50$Variable[
    which(quantile.50$Variable.importance >= 0.1)]) # seuil à modifier 
}

# En plus de l'importance des variables, il faut aussi regarder la forme des
# courbes de réponse pour éliminer les variables / modèles qui donnent des
# réponses qui n'ont pas de sens d'un point de vue biologique

saveRDS(sel_vars, file = paste0("data/selected_variables.RDS"))

