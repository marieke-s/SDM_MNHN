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
# Setting up

# Load libraries
library(biomod2) # Testé fonctionnel en v. 4.2-5-2 --> pb de rétrocompatibilité avec les versions antérieures possible
library(ggtext)
library(RColorBrewer) # Pour les palettes de couleurs
library(terra)
library(tidyterra)

# Style du script
styler::style_file(rstudioapi::getActiveDocumentContext()$path)

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
    proj_names$year_lu,
    sep = "_"
  )
)

# --------------> Important : variables sélectionnées par espèce
#
# Chaque espèce peut avoir ses propres variables sélectionnées. En
# programmation, le plus simple est d'avoir les variables stockées sous forme
# de liste.
# Par exemple, on peut créer la liste des variables sélectionnées comme ceci :
sel_vars <- list(
  Jaime.lannisterii = c("bio5", "bio6", "forests"),
  Tyrion.lannisterii = c("bio5", "bio6", "forests"),
  Jon.snowii = c("bio5", "bio6", "forests"),
  Ned.starkii = c("bio5", "bio6", "forests")
)

# On va utiliser ce format dans ce script, à vous d'adapter la liste en
# adoptant ce format pour vos espèces et vos variables.

# Si vous utilisez le script de sélection de variables (script 7), dans ce cas
# la liste des variables sélectionnées par espèce sera automatiquement
# stockée sur le disque, et on peut alors la charger avec cette ligne de
# commande :

# sel_vars <- readRDS("data/selected_variables.RDS")



# On va faire une boucle qui va lancer toutes les étapes de la modélisation
# pour toutes les espèces, une par une
for (i in 1:nrow(sp_list))
{
  # Pour tester la boucle:
  #i <- 3

  # Nom de l'espèce
  sp <- sp_list$sp[i]

  # Points d'occurrence de l'espèce # Si on travail avec un df plutot que raster, P_points = spp + pred
  P_points <- readRDS(paste0("data/occurrences_", sp))

  # Variables sélectionnes pour l'espèce
  cur_vars <- sel_vars[[sp]]

  # 1. Biomod2 doit-il générer des pseudoabsences/background ?
  if (sp_list$pa.generation[i] == "yes") {
    nb_PA <- 1000 # Normalement il en faudrait beaucoup plus, mais plus on en met plus le model est long à tourner
    runs_PA <- 2 # Nombre de répétitions de pseudo-absences # car c'est de l'aléatoire donc il faut une répétition, si on avait généré 30000 pts la répétition n'aurait pas été nécessaire.
    PA_strategy <- "random" # Stratégie de génération des pseudo-absences
  } else {
    runs_PA <- 0
    nb_PA <- 0
    PA_strategy <- NULL
  }

  # 2. Les données environnementales viennent d'un raster ou d'une matrice ?
  # Dans votre cas actuel tout est stocké en objet spatial (i.e. rasters)
  if (sp_list$env.data.type[i] == "Spatial") {
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

  # 4. Préparation d'un sous-dossier pour notre espèce
  if (!dir.exists(paste0("models/", sp))) # Si le dossier n'existe pas, on le crée
    {
      dir.create(paste0("models/", sp), recursive = T)
    }

  # 5. Fonction d'initialisation de biomod = formatage des données pour biomod avec BIOMOD_FormatingData()
  run_data <- BIOMOD_FormatingData(
    resp.name = sp, # Nom de l'espèce
    resp.var = P_points, # Variable réponse (occurrences)
    expl.var = calib_env_data, # Variables explicatives
    dir.name = "models", # Dossier dans lequel les modèles seront écrits
    resp.xy = coorxy, # Coordonnées d'occurrence
    # Il est possible d'ajouter un jeu de données d'évaluation ici, avec les
    # arguments eval.XXXX
    PA.nb.rep = runs_PA, # Nombre de runs de pseudoabsence
    PA.nb.absences = nb_PA, # Nombre de pseudoabsences
    PA.strategy = PA_strategy
  ) # Stratégie de sélection de pseudoabsences

  # Sauvergarde (dans un dossier que je connais - car biomod save dans des dossiers obscures)
  saveRDS(run_data, file = paste0("models/", sp, "/run_data.RDS"))

  # Check run_data
  run_data
  str(run_data, max.level = 3) # max.level = 3 pour ne pas afficher les listes

  # Illustration des pseudoabsences générées par biomod
  plot(run_data)

  # 6. Calibration des modèles

  # Etapes préalables à envisager :
  # - Cross-validation (block CV, cf. script 8; ou voir BIOMOD_CrossValidation())
  # - Model tuning: explorer la fonction BIOMOD_tuning()
  #
  # Paramétrisation du random forest
  # rf_size <- round(0.7 * length(P_points))
  #
  # a <- BIOMOD_ModelingOptions(
  #   RF = list(ntree = 1000,
  #             sampsize = c("0" = rf_size, "1" = rf_size),
  #             replace = TRUE))


  model_runs <- BIOMOD_Modeling(
    bm.format = run_data, # Données initialisées par biomod
    modeling.id = "1", # A enlever si vous voulez écrire chaque run dans des
    # dossiers différents
    models = c(
      "GLM", # Liste des modèles à lancer
      "GAM",
      "ANN",
      "MARS",
      "GBM",
      "FDA",
      "RF", # cf. Valavi 2021 pour paramétrisation
      "MAXNET",
      "XGBOOST"
    ),
    CV.strategy = "random", # Stratégie de validation croisée
    CV.nb.rep = 2, # Nombre de runs de validation croisée
    CV.perc = 0.7, # ratio gardé pour la calibration dans la CV
    CV.user.table = NULL, # Table de validation croisée (utile si block CV)
    CV.do.full.models = FALSE, # Faut-il faire des modèles finaux avec 100% des
    # données ?
    OPT.data.type = "binary", # Type de données en input
    OPT.strategy = "bigboss", # Quels paramètres pour les modèles ? Les paramètre bigboss sont les param par défaut choisi par les développeurs de biomod, si "default" c'est les param par défaut de chaque modèle (donc chaque package). --> le mieux c'est de faire "user.defined" et de définir les paramètres soi-même.
    # on peut voir les paramètres par défaut dans help documentation.
    # default, bigboss, tuned ou user.defined
    OPT.user.val = NULL, # Vos paramètres customisés pour les modèles ici
    OPT.user.base = "default", # Quels paramètres par défaut pour les modèles
    # pour les paramètres que vous n'avez pas ajustés ?
    weights = NULL, # Entrer manuellement les poids des observations
    # weights ou prevalence au choix # a utiliser que si on prepare toutes les donnees PA nous même
    prevalence = 0.5, # Prévalence des poids entre présence et pseudo-absences
    # ou backgrounds. Notez que ce paramètre ne fonctionne pas correctement
    # pour tous les modèles (notamment le RF) # réduit l'overfitting quand desequilibre de classe
    var.import = 3, # Nombre de runs de randomisation de variable importance
    metric.eval = c("TSS", "ROC", "BOYCE"), # Metriques d'évaluation --> pour présence seule on ne peut utiliser que du Boyce
    nb.cpu = 4, # Pour paralléliser
    do.progress = TRUE
  ) # Faut-il afficher la barre de progression des calculs ?

  # Sauvegarde des modèles calibrés sur le disque
  saveRDS(model_runs, file = paste0("models/", sp, "/model_runs.RDS"))

  # Exploration du contenu
  str(model_runs, max.level = 3)

  # Accéder au contenu : 2 méthodes pour accéder à l'importance des variables :
  # Avec la fonctions "get_...()"
  get_variables_importance(model_runs)

  #  Avec get(load()) sur les liens dans l'objet
  get(load(model_runs@variables.importance@link)) # L'objet est aussi chargé en
  # mémoire dans R (ici "objValue")

  # 6. Préparation du modèle d'ensemble

  # Le modèle d'ensemble tel qu'effectué par biomod présente des limites.
  # Il vaut mieux le construire soi-même
  # plutôt que d'utiliser biomod, par exemple pour évaluer les incertitudes
  # en profondeur ou pour se passer des limitations de cette fonction
  # (i.e. ne permet pas d'inclure des répétitions hors biomod telles que des
  # pseudoabsences générées manuellement,
  # nombre de métriques limitées, ne calcule pas les écarts-types...)
  em_runs <- BIOMOD_EnsembleModeling(
    bm.mod = model_runs, # Modèles individuels calibrés
    models.chosen = "all", # Filtrage des modèles pour l'EM
    em.by = "all", # Sur quelle base les modeles choisis sont combinés (all, algo, PA, PA+algo, PA+run), ici all = on inclus tous les modèles et toutes leurs versions.
    em.algo = c("EMmean", "EMmedian", "EMcv", "EMci", "EMca", "EMwmean"), # Quelles méthodes pour le modèle d'ensemble ? # EMwmean moyenne pondérée par l'évaluation = non conseillé car sous-tend qu'on a suffisamment confiance dans notre métrique d'évaluation. # EMca : somme par pixel et on décide si selon une règle de démocration on juge présence ou abscence (résultat binaire 0 ou 1), EMci = donne l'intervalle de confiance (sort 2 carte ci Inf ate ci Max, pour les 2 bornes).
    metric.select = "BOYCE", # Quelle métrique utiliser pour filtrer les mauvais modèles ou pondérer la contribution des modèles dans l'EM ?
    metric.select.thresh = 0.5, # Quel seuil de filtration des mauvais modèles ? (Boyce : plus il est proche de 1 plus il est informatif, seuil : pas de règle prédéfini, ne pas descendre en dessous de 0.3 = proche de l'aléatoire). Probleme avec Boyce : nécessite au moins 20/30 occurences pour qu'il fasse du sens.
    metric.eval = c("TSS", "ROC", "BOYCE"), # Quelles métriques utiliser pour évaluer l'EM ? --> Si vous n'avez pas de jeu de donné indep vous ne pouvez pas l'évaluer --> l'évaluatoion vient de l'évaluation des modèles individuels. --> evaluation du modèle d'ensemble avec les courbes de réponses et les cartes --> = évaluer la vraisemblence écologique. --> on peut quand même faire l'évaluation sur des données non indep. pour vérifier que l'évaluation n'est pas mauvais : si mauvaise c'est que qqchose ne va vraiment pas.
    var.import = 1, # Runs de randomisation des variables pour calculer leur importance pour l'EM
    prob.ci.alpha = 0.05, # Seuil de l'intervalle de confiance
    nb.cpu = 4
  ) # Pour paralléliser

  # Sauvegarder l'EM sur le disque
  saveRDS(em_runs, file = paste0("models/", sp, "/em_runs.RDS"))

  # 7. Projection des modèles individuels
  # Pour présent = 36 cartes produites

  for (j in list_stacks)
  {
    # Heure du début de la projection
    cat(paste("---- ", Sys.time(), "Projection:", j, "----"))

    # On charge les données environnementales sur lesquelles on projette
    projection_stack <- rast(paste0("data/", j, ".tif"))

    # On projette les modèles
    projection_runs <- BIOMOD_Projection(
      bm.mod = model_runs, # Modèles calibrés
      proj.name = j, # Nom de la projection de l'itération actuelle
      new.env = projection_stack[[cur_vars]], # Données environnementales sur
      # lesquelles on projette les modèles
      models.chosen = "all", # Modèles à projeter
      metric.binary = "TSS", # Avec quelle métrique faut-il transformer les
      # probas en présence-absence ? # Attention TSS n'est pas la meilleure méthode pour P-A et à ne pas utiliser pour les données de présence seule. # Malgrès tout assez répendu dand le domaine.
      metric.filter = NULL, # Avec quelle métrique appliquer un seuil en dessous
      # duquel la proba est forcée à zéro ? # jamais utilisée --> pas conseillé !
      build.clamping.mask = TRUE, # Le clamping mask illustre les zones où les
      # prédictions sont en dehors des valeurs utilisées lors de la calibration
      # Par exemple quand on projette dans le futur et qu'on rencontre des valeurs jamais rencontrées dans la calibration --> permet d'identifier où il y'a de l'extrapolation.
      nb.cpu = 4 # Parallélisation
    )

    saveRDS(projection_runs, file = paste("models/", sp, "/proj_", j,
      "/projection_runs.RDS",
      sep = ""
    ))

    # 8. Projection de l'EM
    proj_em <- BIOMOD_EnsembleForecasting(
      bm.em = em_runs, # Objet correspond à l'EM de biomod
      bm.proj = projection_runs, # Objet issu des projections des modèles
      # individuels
      metric.binary = "TSS"
    ) # Avec quelle métrique faut-il transformer les
    # probas en présence-absence ?

    saveRDS(proj_em, file = paste("models/", sp, "/proj_", j,
      "/projection_em.RDS",
      sep = ""
    ))
    cat(paste("---- ", Sys.time(), "Projection:", j, "finished ----\n"))
    # print(warnings())

    # Vidage de la mémoire
    # rm(list = c("projection_runs", "proj_em"))
  }
}
