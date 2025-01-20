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


library(biomod2) # Testé fonctionnel en version 4.2-5-2
library(RColorBrewer)
library(terra)

# Données baseline 
baseline <- rast("data/baseline.tif")

# Liste des espèces
# Ce tableau contient les paramètres de modélisation :
# nature des données d'occurrence (presence absence ou presence only)
# nature des données environnementales (raster ou data.frame)
# nécessité de générer des points de background / pseudoabsences ou non													
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")

# Tableau avec les noms des projections
proj_names <- readRDS("data/projection_names.RDS")

# On crée un vecteur avec les noms finaux de nos projections
list_stacks <- c(
  "baseline", 
  paste(proj_names$gcm_lu,
        proj_names$scenar_clim,
        proj_names$year_lu, sep = "_")
)


# --------------> Important : variables sélectionnées par espèce
#  
# Chaque espèce peut avoir ses propres variables sélectionnées. En
# programmation, le plus simple est d'avoir les variables stockées sous forme 
# de liste.
# Par exemple, on peut créer la liste des variables sélectionnées comme ceci :
# sel_vars <- list(
#   Jaime.lannisterii = c("bio5", "bio6", "forests"),
#   Tyrion.lannisterii = c("bio5", "bio6", "forests"),
#   Jon.snowii = c("bio5", "bio6", "forests"),
#   Ned.starkii = c("bio5", "bio6", "forests")
# )
# On va utiliser ce format dans ce script, à vous d'adapter la liste en 
# adoptant ce format pour vos espèces et vos variables.

# Si vous utilisez le script de sélection de variables (script 7), dans ce cas
# la liste des variables sélectionnées par espèce sera automatiquement 
# stockée sur le disque, et on peut alors la charger avec cette ligne de 
# commande :

sel_vars <- readRDS("data/selected_variables.RDS")



# On va faire une boucle qui va lancer toutes les étapes de la modélisation
# pour toutes les espèces, une par une
for (i in 1:nrow(sp_list))
{
  # Pour tester la boucle:
  # i <- 3
  
  # Nom de l'espèce
  sp <- sp_list$sp[i]


  # Variables sélectionnes pour l'espèce
  cur_vars <- sel_vars[[sp]]
  
  # 1. Y a-t-il des pseudoabsences/background ?
  if(sp_list$pa.generation[i] == "yes")
  {
    PA_table <- readRDS(paste0("data/bg_table_", sp))
    PA_strategy <- "user.defined"
    # Points d'occurrence de l'espèce
    P_points <- readRDS(paste0("data/occurrences_bg_", sp))  
  } else
  {
    PA_table <- NULL
    PA_strategy <- NULL
    # Points d'occurrence de l'espèce
    P_points <- readRDS(paste0("data/occurrences_", sp))  
  }
  
  # 2. Les données environnementales viennent d'un raster ou d'une matrice ?
  # Dans votre cas actuel tout est stocké en objet spatial (i.e. rasters)
  if(sp_list$env.data.type[i] == "Spatial")
  { 
    calib_env_data <- baseline[[cur_vars]]
  } else
  # Cas où les données environnementales sont stockées sous forme de matrice
  {
    calib_env_data <- P_points[, cur_vars]
  }
  
  # 3. Formatage des occurrences pour biomod
  # On crée un objet avec les coordonnées XY
  coorxy <- P_points[, c("x", "y")]
  # Et un objet avec l'information présence / absence
  P_points <- P_points[, "Observed"]
  
  # 3b. Chargement des folds de cross-validation
  folds <- readRDS(paste0("./data/folds_", sp, ".RDS"))
  
  
  # 4. Préparation d'un sous-dossier pour notre espèce
  if(!dir.exists(paste0("models/", sp)))
  {
    dir.create(paste0("models/", sp), recursive = T)
  }
  
  
  
  # 5. Fonction d'initialisation de biomod
  run_data <- BIOMOD_FormatingData(
    resp.name = sp, # Nom de l'espèce
    resp.var = P_points, # Variable réponse (occurrences)
    expl.var = calib_env_data, # Variables explicatives
    dir.name = "models", # Dossier dans lequel les modèles seront écrits
    resp.xy = coorxy, # Coordonnées d'occurrence
    # Il est possible d'ajouter un jeu de données d'évaluation ici, avec les
    # arguments eval.XXXX
    PA.nb.rep = ncol(PA_table),
    PA.strategy = PA_strategy, 
    PA.user.table = PA_table)
  
  saveRDS(run_data, file = paste0("models/", sp, "/run_data.RDS"))
  
  # Exploration du contenu
  run_data
  str(run_data, max.level = 3)
  
  
  # Illustration des pseudoabsences générées par biomod
  plot(run_data)

  # 6. Calibration des modèles
  
  ##### Paramétrisation des modèles #####
  # La paramétrisation des modèles est importante, mais elle est très complexe
  # à mettre en place dans biomod2 actuellement (01/2024, v4.2-5),
  # pour deux raisons :
  #  - la documentation n'est pas à jour / disponible
  #  - il faut bien respecter les conventions de nommage des objets
  #  - chaque run (cross-validation + background) doit être paramétré 
  #    individuellement
  #  - les messages d'erreurs ou d'information n'indiquent pas clairement si votre
  #    paramétrisation a été prise en compte
  # 
  # Il faut donc être très attentif à ce que vous faites et aux résultats des 
  # modèles.
  # Voici donc une proposition fonctionnelle sur la v 4.2-5 de biomod2, mais
  # nécessairement complexe en termes de programmation.
  
  # Nous allons modifier deux types de paramètres ici :
  # 1. les poids des présences et absences/background. Nous allons manuellement
  # assigner des poids de telle sorte que la somme du poids des présences soit
  # égale à la somme du poids des absences
  # 2. appliquer du down-sampling sur le random-forest, afin de s'assurer qu'il
  # y a le même nombre de présences et absences/background dans chaque arbre
  # individuel calibré par le RF. (cf. Valavi et al. 2021 Ecography)
  # Ainsi que quelques autres paramètres en s'inspirant des méthodes utilisées
  # par Valavi et al. 2021 (Ecological Monographs)
  #
  # Ces deux exemples vous donneront un bon point de départ pour commencer,
  # à vous de modifier les paramètres de manière plus complexe et complète si
  # vous le souhaitez : regardez quels paramètres sont modifiables pour chaque
  # modèle en inspectant l'aide des packages de chaque modèle implémenté
  # dans biomod2. Pour connaître les packages et fonctions de modélisation,
  # tapez ModelsTable dans la console.
  
  
  # Tout d'abord, on va récupérer le nombre de présences et absences/bg
  # par run
  calib_summary <- summary(run_data, calib.lines =  folds) 
  # Bien entendu, on ne travaille que sur les données de calibration ici, 
  # car il s'agit de paramétrer la calibration des modèles :
  calib_summary <- calib_summary[
    which(calib_summary$data == "calibration"), ]
  
  
  # Il va falloir paramétrer les modèles pour chacun des runs, soit chaque 
  # ligne du tableau calib_summary
  
  # Nous allons conserver un objet par modèle pour les paramètres.
  # Ajoutez/modifiez cette liste en fonction des modèles que vous prévoyez
  # de faire
  GLM_param_list <- NULL
  GAM_param_list <- NULL
  GBM_param_list <- NULL
  MARS_param_list <- NULL
  RF_param_list <- NULL
  MAXNET_param_list <- NULL
  XGBOOST_param_list <- NULL
  
  for (cvrun in 1:nrow(calib_summary)) {
    
    ### 1. On récupère les infos sur le run actuel : occurrences, nombre
    # de présences, nombre de backgrounds/absences
    
    occurrences <- run_data@data.species[
      # On récupère les occurrences pour le run actuel :
      # PA et CV
      which(folds[, paste0("_", calib_summary$PA[cvrun],
                   "_", calib_summary$run[cvrun])])
    ]
    
    # On va récupérer le nombre de présences pour la calibration ici
    prNum <- calib_summary$Presences[cvrun]
    # Et le nombre d'absences ou bg ici
    if(sp_list$pa.generation[i] == "no")
    {
      bgNum <- calib_summary$True_Absences[cvrun]
    } else {
      bgNum <- calib_summary$Pseudo_Absences[cvrun]
    }
    
    # Nous allons attribuer un poids aux occurrences, de telle sorte que s'il 
    # s'agit d'une présence, le poids est 1, et si c'est une absence / un point
    # de background, le poids est égal au ratio présence / absence-bg
    # Cette formule permet de s'assurer que la somme du poids des présences est
    # égale à la somme du poids des absences
    
    wt <- ifelse(occurrences == 1, # Si c'est une présence 
                 1, # On met un poids de 1
                 prNum / bgNum) # Sinon on met un poids égal au ratio pres/abs
    
    
    ### 2. On complète les listes de paramètres pour chaque modèle avec les
    # infos sur le run actuel
    
    # GLM :
    # On va corriger le déséquilibre de classe en attribuant des poids aux
    # présences et aux absences/background 
    # Pour paramétrer un GLM finement il faut passer par une procédure
    # stepwise, le mieux étant de fitter une régression régularisée,
    # mais cette méthode n'est pas disponible dans biomod2
    # cf. https://rvalavi.github.io/Presence-Only-SDM/#glm et section suivante
    GLM_param_list[[paste0("_",
                           calib_summary$PA[[cvrun]],
                           "_",
                           calib_summary$run[[cvrun]])]] =
      list(weights = wt)
    
    
    # GBM :
    # On va corriger le déséquilibre de classe en attribuant des poids aux
    # présences et aux absences/background. On va utiliser un grand nombre 
    # d'arbres (500), et on va utiliser les paramètres par défaut dérivés
    # des tests de Valavi et al. (2021, Ecological Monographs)
    # https://rvalavi.github.io/Presence-Only-SDM/#brt-gbm
    GBM_param_list[[paste0("_",
                           calib_summary$PA[[cvrun]],
                           "_",
                           calib_summary$run[[cvrun]])]] =
      list(interaction.depth = 5,
           n.trees = 500, 
           shrinkage = 0.001,
           bag.fraction = 0.75,
           cv.folds = 5,
           weights = wt)
    
    # GAM :
    # On va corriger le déséquilibre de classe en attribuant des poids aux
    # présences et aux absences/background. On les paramètres par défaut dérivés
    # des tests de Valavi et al. (2021, Ecological Monographs)
    # https://rvalavi.github.io/Presence-Only-SDM/#gam
    GAM_param_list[[paste0("_",
                           calib_summary$PA[[cvrun]],
                           "_",
                           calib_summary$run[[cvrun]])]] <-     
      list(weights = wt,
           k = 10,
           method = "REML")
    
    # MARS :
    # On va corriger le déséquilibre de classe en attribuant des poids aux
    # présences et aux absences/background. Pour ajuster les autres paramètres
    # du modèle, il faut passer par une procédure de validation croisée,
    # cf. https://rvalavi.github.io/Presence-Only-SDM/#mars
    MARS_param_list[[paste0("_",
                            calib_summary$PA[[cvrun]],
                            "_",
                            calib_summary$run[[cvrun]])]] <- 
      list(weights = wt)
  
    
    # Random forest :
    # On fait du down-sampling, i.e. s'assurer qu'il y a autant de présences
    # que d'absences dans les arbres individuels calibrés en interne par le
    # random forest. On va aussi permettre de rééchantillonner les mêmes points
    # pour un même arbre. On met un nombre d'arbres assez grand (1000 ici). Ces
    # paramètres sont tirés de Valavi et al. 2021, Ecography et Ecological 
    # Monographs
    # https://rvalavi.github.io/Presence-Only-SDM/#rf-and-rf-down-sampled
    RF_param_list[[paste0("_",
                          calib_summary$PA[[cvrun]],
                          "_",
                          calib_summary$run[[cvrun]])]] =
      list(ntree = 1000,
           sampsize =  c("0" = prNum,
                         "1" = prNum),
           replace = TRUE)
    

    # MAXNET
    # Les valeurs par défaut de MaxNet semble donner de très bons résultats 
    # donc nous gardons ces paramètres par défaut ici, en ajustant simplement
    # le paramètre de régularisation à 1
    # https://rvalavi.github.io/Presence-Only-SDM/#maxent-and-maxnet
    MAXNET_param_list[[paste0("_",
                            calib_summary$PA[[cvrun]],
                            "_",
                            calib_summary$run[[cvrun]])]] <- 
      list(regmult = 1)
    
    # XGBOOST
    # On va utiliser des paramètres similaires à Valavi et al. 2021 pour 
    # XGBOOST, sachant que ces paramètres n'ont pas donné de très bons 
    # résultats. Il faudrait donc affiner la paramétrisation de ce modèle.
    # https://rvalavi.github.io/Presence-Only-SDM/#xgboost
    XGBOOST_param_list[[paste0("_",
                               calib_summary$PA[[cvrun]],
                               "_",
                               calib_summary$run[[cvrun]])]] <-
      list(nrounds = 5000,
           eta = 0.001,
           max_depth = 5,
           subsample = 0.75,
           gamma = 0,
           colsample_bytree = 0.8,
           min_child_weight = 1,
           weight = wt,
           verbose = 0)
  }
  
  
  model_parameters <- bm_ModelingOptions(
    data.type = "binary", # Données d'occurrence = binary
    models = c("GLM", 
               "GAM.mgcv.gam", # biomod2 a plusieurs GAM donc il faut préciser
               # Cf. ModelsTable
               "GBM", 
               "MARS",
               "RF", 
               "MAXNET", 
               "XGBOOST"),
    strategy = "user.defined", # On va définir les paramètres nous-mêmes
    user.base = "bigboss", # Quels paramètres par défaut biomod doit-il utiliser
    # pour tous les paramètres que l'on n'a pas customisés ? 
    user.val = list(
      GLM.binary.stats.glm = GLM_param_list,
      GBM.binary.gbm.gbm = GBM_param_list,
      GAM.binary.mgcv.gam = GAM_param_list,
      MARS.binary.earth.earth = MARS_param_list,
      RF.binary.randomForest.randomForest = RF_param_list,
      XGBOOST.binary.xgboost.xgboost = XGBOOST_param_list
    ),
    bm.format = run_data,
    calib.lines = folds
  )
  # biomod2 va vous mettre de nombreux messages d'avis en v4.2-5, n'en tenez pas
  # compte, ils sont inappropriés
  # cf. https://github.com/biomodhub/biomod2/issues/402
  
  
  
  
  
  model_runs <- BIOMOD_Modeling(
    bm.format = run_data, # Données initialisées par biomod
    modeling.id = "1", # A enlever si vous voulez écrire chaque run dans des 
    # dossiers différents
    models =  c("GLM", # Liste des modèles à lancer
                "GAM",
                "MARS", 
                "GBM",
                "RF",
                "MAXNET",
                "XGBOOST"),
    CV.strategy = "user.defined", # Stratégie de validation croisée
    # CV.nb.rep = 2, # Nombre de runs de validation croisée
    # CV.perc = 0.7, # ratio gardé pour la calibration dans la CV
    CV.user.table = folds, # Table de validation croisée (utile si block CV)
    CV.do.full.models = FALSE, # Faut-il faire des modèles finaux avec 100% des
    # données ?
    OPT.data.type = "binary", # Type de données en input
    OPT.strategy = "default", # Quels paramètres pour les modèles ?
    OPT.user.val = model_parameters, # Vos paramètres customisés 
    OPT.user.base = "default", # Quels paramètres par défaut pour les modèles
    # pour les paramètres que vous n'avez pas ajustés ?
    weights = NULL, # Entrer manuellement les poids des observations
    prevalence = 0.5, # Prévalence des poids entre présence et pseudo-absences
    # ou backgrounds. Notez que ce paramètre ne fonctionne pas correctement
    # pour tous les modèles (notamment le RF - mais on l'a paramétrisé 
	# manuellement pour corriger ce problème)
    var.import = 3, # Nombre de runs de randomisation de variable importance
    metric.eval = c("TSS", "ROC", "BOYCE"), # Metriques d'évaluation
    nb.cpu = 4, # Pour paralléliser 
    do.progress = TRUE) # Faut-il afficher la barre de progression des calculs ?
  
  # Sauvegarde des modèles calibrés sur le disque
  saveRDS(model_runs, file = paste0("models/", sp, "/model_runs.RDS"))

  # Exploration du contenu
  str(model_runs, max.level = 3)
  # Accéder au contenu : deux méthodes
  # Exemple de l'importance des variables
  # Fonctions "get_...()"
  get_variables_importance(model_runs)
  # Charger avec les liens dans l'objet
  get(load(model_runs@variables.importance@link)) # L'objet est aussi chargé en mémoire dans R (ici "objValue")

  # 6. Préparation du modèle d'ensemble
  
  # Le modèle d'ensemble tel qu'effectué par biomod 
  # Il vaut mieux le construire soi-même
  # plutôt que d'utiliser biomod, par exemple pour évaluer les incertitudes
  # en profondeur ou pour se passer des limitations de cette fonction 
  # (i.e. ne permet pas d'inclure des répétitions hors biomod telles que des
  # pseudoabsences générées manuellement,
  # nombre de métriques limitées, ne calcule pas les écarts-types...)
  em_runs <- BIOMOD_EnsembleModeling(
    bm.mod = model_runs, # Modèles individuels calibrés
    models.chosen = 'all', # Filtrage des modèles pour l'EM
    em.by = 'all', # Qu'est-ce qui doit être inclus dans l'EM
    em.algo = c("EMmean", "EMmedian", "EMcv", "EMci", "EMca", "EMwmean"), #
    # Quelles méthodes pour le modèle d'ensemble ?
    metric.select = 'BOYCE', # Quelle métrique utiliser pour filtrer les
    # mauvais modèles 
    # ou pondérer la contribution des modèles dans l'EM ?
    metric.select.thresh = 0.5, # Quel seuil de filtration des mauvais modèles ? 
    metric.eval = c("TSS", "ROC", "BOYCE"), # Quelles métriques utiliser pour é
    # valuer l'EM ?
    var.import = 1, # Runs de variable importance pour l'EM
    prob.ci.alpha = 0.05, # Seuil de l'intervalle de confiance
    nb.cpu = 4) # Pour paralléliser 
  
  # Sauvegarder l'EM sur le disque
  saveRDS(em_runs, file = paste0("models/", sp, "/em_runs.RDS"))
c
  for(j in list_stacks)
  {
    # 7. Projection des modèles individuels
    
    cat(paste("---- ", Sys.time(), "Projection:", j, "----"))
    projection_stack <- rast(paste0("data/", j, ".tif"))
    # Dans ce cas, biomod n'aimera pas et ne fera pas la projection
    projection_runs <- BIOMOD_Projection(
      bm.mod = model_runs, # Modèles calibrés
      proj.name = j, # Nom de la projection actuelle
      new.env = projection_stack[[cur_vars]], # Données environnementales sur
      # lesquelles on projette les modèles
      models.chosen = "all", # Modèles à projeter
      metric.binary = "TSS", # Avec quelle métrique faut-il transformer les 
      # probas en présence-absence ? 
      metric.filter = NULL, # Avec quelle métrique appliquer un seuil en dessous
      # duquel la proba est forcée à zéro ?
      build.clamping.mask = TRUE, # Le clamping mask illustre les zones où les 
      # prédictions sont en dehors des valeurs
      # utilisées lors de la calibration
      nb.cpu = 4) # Parallélisation
    saveRDS(projection_runs, file = paste("models/", sp, "/proj_", j, 
                                          "/projection_runs.RDS", sep = ""))
    
    # 8. Projection de l'EM
    proj_em <- BIOMOD_EnsembleForecasting(
      bm.em = em_runs, # Objet correspond à l'EM de biomod
      bm.proj = projection_runs, # Objet issu des projections des modèles
      # individuels
      metric.binary = "TSS") # Avec quelle métrique faut-il transformer les 
    # probas en présence-absence ? 
    
    saveRDS(proj_em, file = paste("models/", sp, "/proj_", j, 
                                  "/projection_em.RDS", sep = ""))
    cat(paste("---- ", Sys.time(), "Projection:", j, "finished ----\n"))
    # print(warnings())
    
    # Vidage de la mémoire
    # rm(list = c("projection_runs", "proj_em"))
  }
}


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


library(biomod2) # Testé fonctionnel en version 4.2-5-2
library(RColorBrewer)
library(terra)

# Données baseline 
baseline <- rast("data/baseline.tif")

# Liste des espèces
# Ce tableau contient les paramètres de modélisation :
# nature des données d'occurrence (presence absence ou presence only)
# nature des données environnementales (raster ou data.frame)
# nécessité de générer des points de background / pseudoabsences ou non													
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")

# Tableau avec les noms des projections
proj_names <- readRDS("data/projection_names.RDS")

# On crée un vecteur avec les noms finaux de nos projections
list_stacks <- c(
  "baseline", 
  paste(proj_names$gcm_lu,
        proj_names$scenar_clim,
        proj_names$year_lu, sep = "_")
)


# --------------> Important : variables sélectionnées par espèce
#  
# Chaque espèce peut avoir ses propres variables sélectionnées. En
# programmation, le plus simple est d'avoir les variables stockées sous forme 
# de liste.
# Par exemple, on peut créer la liste des variables sélectionnées comme ceci :
# sel_vars <- list(
#   Jaime.lannisterii = c("bio5", "bio6", "forests"),
#   Tyrion.lannisterii = c("bio5", "bio6", "forests"),
#   Jon.snowii = c("bio5", "bio6", "forests"),
#   Ned.starkii = c("bio5", "bio6", "forests")
# )
# On va utiliser ce format dans ce script, à vous d'adapter la liste en 
# adoptant ce format pour vos espèces et vos variables.

# Si vous utilisez le script de sélection de variables (script 7), dans ce cas
# la liste des variables sélectionnées par espèce sera automatiquement 
# stockée sur le disque, et on peut alors la charger avec cette ligne de 
# commande :

sel_vars <- readRDS("data/selected_variables.RDS")



# On va faire une boucle qui va lancer toutes les étapes de la modélisation
# pour toutes les espèces, une par une
for (i in 1:nrow(sp_list))
{
  # Pour tester la boucle:
  # i <- 3
  
  # Nom de l'espèce
  sp <- sp_list$sp[i]


  # Variables sélectionnes pour l'espèce
  cur_vars <- sel_vars[[sp]]
  
  # 1. Y a-t-il des pseudoabsences/background ?
  if(sp_list$pa.generation[i] == "yes")
  {
    PA_table <- readRDS(paste0("data/bg_table_", sp))
    PA_strategy <- "user.defined"
    # Points d'occurrence de l'espèce
    P_points <- readRDS(paste0("data/occurrences_bg_", sp))  
  } else
  {
    PA_table <- NULL
    PA_strategy <- NULL
    # Points d'occurrence de l'espèce
    P_points <- readRDS(paste0("data/occurrences_", sp))  
  }
  
  # 2. Les données environnementales viennent d'un raster ou d'une matrice ?
  # Dans votre cas actuel tout est stocké en objet spatial (i.e. rasters)
  if(sp_list$env.data.type[i] == "Spatial")
  { 
    calib_env_data <- baseline[[cur_vars]]
  } else
  # Cas où les données environnementales sont stockées sous forme de matrice
  {
    calib_env_data <- P_points[, cur_vars]
  }
  
  # 3. Formatage des occurrences pour biomod
  # On crée un objet avec les coordonnées XY
  coorxy <- P_points[, c("x", "y")]
  # Et un objet avec l'information présence / absence
  P_points <- P_points[, "Observed"]
  
  # 3b. Chargement des folds de cross-validation
  folds <- readRDS(paste0("./data/folds_", sp, ".RDS"))
  
  
  # 4. Préparation d'un sous-dossier pour notre espèce
  if(!dir.exists(paste0("models/", sp)))
  {
    dir.create(paste0("models/", sp), recursive = T)
  }
  
  
  
  # 5. Fonction d'initialisation de biomod
  run_data <- BIOMOD_FormatingData(
    resp.name = sp, # Nom de l'espèce
    resp.var = P_points, # Variable réponse (occurrences)
    expl.var = calib_env_data, # Variables explicatives
    dir.name = "models", # Dossier dans lequel les modèles seront écrits
    resp.xy = coorxy, # Coordonnées d'occurrence
    # Il est possible d'ajouter un jeu de données d'évaluation ici, avec les
    # arguments eval.XXXX
    PA.nb.rep = ncol(PA_table),
    PA.strategy = PA_strategy, 
    PA.user.table = PA_table)
  
  saveRDS(run_data, file = paste0("models/", sp, "/run_data.RDS"))
  
  # Exploration du contenu
  run_data
  str(run_data, max.level = 3)
  
  
  # Illustration des pseudoabsences générées par biomod
  plot(run_data)

  # 6. Calibration des modèles
  
  ##### Paramétrisation des modèles #####
  # La paramétrisation des modèles est importante, mais elle est très complexe
  # à mettre en place dans biomod2 actuellement (01/2024, v4.2-5),
  # pour deux raisons :
  #  - la documentation n'est pas à jour / disponible
  #  - il faut bien respecter les conventions de nommage des objets
  #  - chaque run (cross-validation + background) doit être paramétré 
  #    individuellement
  #  - les messages d'erreurs ou d'information n'indiquent pas clairement si votre
  #    paramétrisation a été prise en compte
  # 
  # Il faut donc être très attentif à ce que vous faites et aux résultats des 
  # modèles.
  # Voici donc une proposition fonctionnelle sur la v 4.2-5 de biomod2, mais
  # nécessairement complexe en termes de programmation.
  
  # Nous allons modifier deux types de paramètres ici :
  # 1. les poids des présences et absences/background. Nous allons manuellement
  # assigner des poids de telle sorte que la somme du poids des présences soit
  # égale à la somme du poids des absences
  # 2. appliquer du down-sampling sur le random-forest, afin de s'assurer qu'il
  # y a le même nombre de présences et absences/background dans chaque arbre
  # individuel calibré par le RF. (cf. Valavi et al. 2021 Ecography)
  # Ainsi que quelques autres paramètres en s'inspirant des méthodes utilisées
  # par Valavi et al. 2021 (Ecological Monographs)
  #
  # Ces deux exemples vous donneront un bon point de départ pour commencer,
  # à vous de modifier les paramètres de manière plus complexe et complète si
  # vous le souhaitez : regardez quels paramètres sont modifiables pour chaque
  # modèle en inspectant l'aide des packages de chaque modèle implémenté
  # dans biomod2. Pour connaître les packages et fonctions de modélisation,
  # tapez ModelsTable dans la console.
  
  
  # Tout d'abord, on va récupérer le nombre de présences et absences/bg
  # par run
  calib_summary <- summary(run_data, calib.lines =  folds) 
  # Bien entendu, on ne travaille que sur les données de calibration ici, 
  # car il s'agit de paramétrer la calibration des modèles :
  calib_summary <- calib_summary[
    which(calib_summary$data == "calibration"), ]
  
  
  # Il va falloir paramétrer les modèles pour chacun des runs, soit chaque 
  # ligne du tableau calib_summary
  
  # Nous allons conserver un objet par modèle pour les paramètres.
  # Ajoutez/modifiez cette liste en fonction des modèles que vous prévoyez
  # de faire
  GLM_param_list <- NULL
  GAM_param_list <- NULL
  GBM_param_list <- NULL
  MARS_param_list <- NULL
  RF_param_list <- NULL
  MAXNET_param_list <- NULL
  XGBOOST_param_list <- NULL
  
  for (cvrun in 1:nrow(calib_summary)) {
    
    ### 1. On récupère les infos sur le run actuel : occurrences, nombre
    # de présences, nombre de backgrounds/absences
    
    occurrences <- run_data@data.species[
      # On récupère les occurrences pour le run actuel :
      # PA et CV
      which(folds[, paste0("_", calib_summary$PA[cvrun],
                   "_", calib_summary$run[cvrun])])
    ]
    
    # On va récupérer le nombre de présences pour la calibration ici
    prNum <- calib_summary$Presences[cvrun]
    # Et le nombre d'absences ou bg ici
    if(sp_list$pa.generation[i] == "no")
    {
      bgNum <- calib_summary$True_Absences[cvrun]
    } else {
      bgNum <- calib_summary$Pseudo_Absences[cvrun]
    }
    
    # Nous allons attribuer un poids aux occurrences, de telle sorte que s'il 
    # s'agit d'une présence, le poids est 1, et si c'est une absence / un point
    # de background, le poids est égal au ratio présence / absence-bg
    # Cette formule permet de s'assurer que la somme du poids des présences est
    # égale à la somme du poids des absences
    
    wt <- ifelse(occurrences == 1, # Si c'est une présence 
                 1, # On met un poids de 1
                 prNum / bgNum) # Sinon on met un poids égal au ratio pres/abs
    
    
    ### 2. On complète les listes de paramètres pour chaque modèle avec les
    # infos sur le run actuel
    
    # GLM :
    # On va corriger le déséquilibre de classe en attribuant des poids aux
    # présences et aux absences/background 
    # Pour paramétrer un GLM finement il faut passer par une procédure
    # stepwise, le mieux étant de fitter une régression régularisée,
    # mais cette méthode n'est pas disponible dans biomod2
    # cf. https://rvalavi.github.io/Presence-Only-SDM/#glm et section suivante
    GLM_param_list[[paste0("_",
                           calib_summary$PA[[cvrun]],
                           "_",
                           calib_summary$run[[cvrun]])]] =
      list(weights = wt)
    
    
    # GBM :
    # On va corriger le déséquilibre de classe en attribuant des poids aux
    # présences et aux absences/background. On va utiliser un grand nombre 
    # d'arbres (500), et on va utiliser les paramètres par défaut dérivés
    # des tests de Valavi et al. (2021, Ecological Monographs)
    # https://rvalavi.github.io/Presence-Only-SDM/#brt-gbm
    GBM_param_list[[paste0("_",
                           calib_summary$PA[[cvrun]],
                           "_",
                           calib_summary$run[[cvrun]])]] =
      list(interaction.depth = 5,
           n.trees = 500, 
           shrinkage = 0.001,
           bag.fraction = 0.75,
           cv.folds = 5,
           weights = wt)
    
    # GAM :
    # On va corriger le déséquilibre de classe en attribuant des poids aux
    # présences et aux absences/background. On les paramètres par défaut dérivés
    # des tests de Valavi et al. (2021, Ecological Monographs)
    # https://rvalavi.github.io/Presence-Only-SDM/#gam
    GAM_param_list[[paste0("_",
                           calib_summary$PA[[cvrun]],
                           "_",
                           calib_summary$run[[cvrun]])]] <-     
      list(weights = wt,
           k = 10,
           method = "REML")
    
    # MARS :
    # On va corriger le déséquilibre de classe en attribuant des poids aux
    # présences et aux absences/background. Pour ajuster les autres paramètres
    # du modèle, il faut passer par une procédure de validation croisée,
    # cf. https://rvalavi.github.io/Presence-Only-SDM/#mars
    MARS_param_list[[paste0("_",
                            calib_summary$PA[[cvrun]],
                            "_",
                            calib_summary$run[[cvrun]])]] <- 
      list(weights = wt)
  
    
    # Random forest :
    # On fait du down-sampling, i.e. s'assurer qu'il y a autant de présences
    # que d'absences dans les arbres individuels calibrés en interne par le
    # random forest. On va aussi permettre de rééchantillonner les mêmes points
    # pour un même arbre. On met un nombre d'arbres assez grand (1000 ici). Ces
    # paramètres sont tirés de Valavi et al. 2021, Ecography et Ecological 
    # Monographs
    # https://rvalavi.github.io/Presence-Only-SDM/#rf-and-rf-down-sampled
    RF_param_list[[paste0("_",
                          calib_summary$PA[[cvrun]],
                          "_",
                          calib_summary$run[[cvrun]])]] =
      list(ntree = 1000,
           sampsize =  c("0" = prNum,
                         "1" = prNum),
           replace = TRUE)
    

    # MAXNET
    # Les valeurs par défaut de MaxNet semble donner de très bons résultats 
    # donc nous gardons ces paramètres par défaut ici, en ajustant simplement
    # le paramètre de régularisation à 1
    # https://rvalavi.github.io/Presence-Only-SDM/#maxent-and-maxnet
    MAXNET_param_list[[paste0("_",
                            calib_summary$PA[[cvrun]],
                            "_",
                            calib_summary$run[[cvrun]])]] <- 
      list(regmult = 1)
    
    # XGBOOST
    # On va utiliser des paramètres similaires à Valavi et al. 2021 pour 
    # XGBOOST, sachant que ces paramètres n'ont pas donné de très bons 
    # résultats. Il faudrait donc affiner la paramétrisation de ce modèle.
    # https://rvalavi.github.io/Presence-Only-SDM/#xgboost
    XGBOOST_param_list[[paste0("_",
                               calib_summary$PA[[cvrun]],
                               "_",
                               calib_summary$run[[cvrun]])]] <-
      list(nrounds = 5000,
           eta = 0.001,
           max_depth = 5,
           subsample = 0.75,
           gamma = 0,
           colsample_bytree = 0.8,
           min_child_weight = 1,
           weight = wt,
           verbose = 0)
  }
  
  
  model_parameters <- bm_ModelingOptions(
    data.type = "binary", # Données d'occurrence = binary
    models = c("GLM", 
               "GAM.mgcv.gam", # biomod2 a plusieurs GAM donc il faut préciser
               # Cf. ModelsTable
               "GBM", 
               "MARS",
               "RF", 
               "MAXNET", 
               "XGBOOST"),
    strategy = "user.defined", # On va définir les paramètres nous-mêmes
    user.base = "bigboss", # Quels paramètres par défaut biomod doit-il utiliser
    # pour tous les paramètres que l'on n'a pas customisés ? 
    user.val = list(
      GLM.binary.stats.glm = GLM_param_list,
      GBM.binary.gbm.gbm = GBM_param_list,
      GAM.binary.mgcv.gam = GAM_param_list,
      MARS.binary.earth.earth = MARS_param_list,
      RF.binary.randomForest.randomForest = RF_param_list,
      XGBOOST.binary.xgboost.xgboost = XGBOOST_param_list
    ),
    bm.format = run_data,
    calib.lines = folds
  )
  # biomod2 va vous mettre de nombreux messages d'avis en v4.2-5, n'en tenez pas
  # compte, ils sont inappropriés
  # cf. https://github.com/biomodhub/biomod2/issues/402
  
  
  
  
  
  model_runs <- BIOMOD_Modeling(
    bm.format = run_data, # Données initialisées par biomod
    modeling.id = "1", # A enlever si vous voulez écrire chaque run dans des 
    # dossiers différents
    models =  c("GLM", # Liste des modèles à lancer
                "GAM",
                "MARS", 
                "GBM",
                "RF",
                "MAXNET",
                "XGBOOST"),
    CV.strategy = "user.defined", # Stratégie de validation croisée
    # CV.nb.rep = 2, # Nombre de runs de validation croisée
    # CV.perc = 0.7, # ratio gardé pour la calibration dans la CV
    CV.user.table = folds, # Table de validation croisée (utile si block CV)
    CV.do.full.models = FALSE, # Faut-il faire des modèles finaux avec 100% des
    # données ?
    OPT.data.type = "binary", # Type de données en input
    OPT.strategy = "default", # Quels paramètres pour les modèles ?
    OPT.user.val = model_parameters, # Vos paramètres customisés 
    OPT.user.base = "default", # Quels paramètres par défaut pour les modèles
    # pour les paramètres que vous n'avez pas ajustés ?
    weights = NULL, # Entrer manuellement les poids des observations
    prevalence = 0.5, # Prévalence des poids entre présence et pseudo-absences
    # ou backgrounds. Notez que ce paramètre ne fonctionne pas correctement
    # pour tous les modèles (notamment le RF - mais on l'a paramétrisé 
	# manuellement pour corriger ce problème)
    var.import = 3, # Nombre de runs de randomisation de variable importance
    metric.eval = c("TSS", "ROC", "BOYCE"), # Metriques d'évaluation
    nb.cpu = 4, # Pour paralléliser 
    do.progress = TRUE) # Faut-il afficher la barre de progression des calculs ?
  
  # Sauvegarde des modèles calibrés sur le disque
  saveRDS(model_runs, file = paste0("models/", sp, "/model_runs.RDS"))

  # Exploration du contenu
  str(model_runs, max.level = 3)
  # Accéder au contenu : deux méthodes
  # Exemple de l'importance des variables
  # Fonctions "get_...()"
  get_variables_importance(model_runs)
  # Charger avec les liens dans l'objet
  get(load(model_runs@variables.importance@link)) # L'objet est aussi chargé en mémoire dans R (ici "objValue")

  # 6. Préparation du modèle d'ensemble
  
  # Le modèle d'ensemble tel qu'effectué par biomod 
  # Il vaut mieux le construire soi-même
  # plutôt que d'utiliser biomod, par exemple pour évaluer les incertitudes
  # en profondeur ou pour se passer des limitations de cette fonction 
  # (i.e. ne permet pas d'inclure des répétitions hors biomod telles que des
  # pseudoabsences générées manuellement,
  # nombre de métriques limitées, ne calcule pas les écarts-types...)
  em_runs <- BIOMOD_EnsembleModeling(
    bm.mod = model_runs, # Modèles individuels calibrés
    models.chosen = 'all', # Filtrage des modèles pour l'EM
    em.by = 'all', # Qu'est-ce qui doit être inclus dans l'EM
    em.algo = c("EMmean", "EMmedian", "EMcv", "EMci", "EMca", "EMwmean"), #
    # Quelles méthodes pour le modèle d'ensemble ?
    metric.select = 'BOYCE', # Quelle métrique utiliser pour filtrer les
    # mauvais modèles 
    # ou pondérer la contribution des modèles dans l'EM ?
    metric.select.thresh = 0.5, # Quel seuil de filtration des mauvais modèles ? 
    metric.eval = c("TSS", "ROC", "BOYCE"), # Quelles métriques utiliser pour é
    # valuer l'EM ?
    var.import = 1, # Runs de variable importance pour l'EM
    prob.ci.alpha = 0.05, # Seuil de l'intervalle de confiance
    nb.cpu = 4) # Pour paralléliser 
  
  # Sauvegarder l'EM sur le disque
  saveRDS(em_runs, file = paste0("models/", sp, "/em_runs.RDS"))
c
  for(j in list_stacks)
  {
    # 7. Projection des modèles individuels
    
    cat(paste("---- ", Sys.time(), "Projection:", j, "----"))
    projection_stack <- rast(paste0("data/", j, ".tif"))
    # Dans ce cas, biomod n'aimera pas et ne fera pas la projection
    projection_runs <- BIOMOD_Projection(
      bm.mod = model_runs, # Modèles calibrés
      proj.name = j, # Nom de la projection actuelle
      new.env = projection_stack[[cur_vars]], # Données environnementales sur
      # lesquelles on projette les modèles
      models.chosen = "all", # Modèles à projeter
      metric.binary = "TSS", # Avec quelle métrique faut-il transformer les 
      # probas en présence-absence ? 
      metric.filter = NULL, # Avec quelle métrique appliquer un seuil en dessous
      # duquel la proba est forcée à zéro ?
      build.clamping.mask = TRUE, # Le clamping mask illustre les zones où les 
      # prédictions sont en dehors des valeurs
      # utilisées lors de la calibration
      nb.cpu = 4) # Parallélisation
    saveRDS(projection_runs, file = paste("models/", sp, "/proj_", j, 
                                          "/projection_runs.RDS", sep = ""))
    
    # 8. Projection de l'EM
    proj_em <- BIOMOD_EnsembleForecasting(
      bm.em = em_runs, # Objet correspond à l'EM de biomod
      bm.proj = projection_runs, # Objet issu des projections des modèles
      # individuels
      metric.binary = "TSS") # Avec quelle métrique faut-il transformer les 
    # probas en présence-absence ? 
    
    saveRDS(proj_em, file = paste("models/", sp, "/proj_", j, 
                                  "/projection_em.RDS", sep = ""))
    cat(paste("---- ", Sys.time(), "Projection:", j, "finished ----\n"))
    # print(warnings())
    
    # Vidage de la mémoire
    # rm(list = c("projection_runs", "proj_em"))
  }
}


>>>>>>> 5768564 (2nd commit)
