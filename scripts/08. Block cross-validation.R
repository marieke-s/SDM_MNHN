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


library(blockCV)
library(terra)

sp_list <- read.csv("./data_cours/species_list.csv", sep = ";")


##### Partie 1 : Préparation des données #####
# 
# La validation croisée par bloc nécessite d'évaluer le range d'autocorrélation
# des variables environnementales. Pour cela, il faut projeter les raster,
# sinnon le calcul de l'autocorrélation spatiale est imprécis.
# Si vos données ne sont pas projetées (e.g. elles sont en WGS84), il faudra
# les projeter au moins pour le calcul du range d'autocorrélation.
# 
# Si vous travaillez sur de larges étendues spatiales couvrant les hautes 
# latitudes, il vaut mieux projeter toutes vos données à l'étape de 
# préparation des données, afin d'obtenir des blocs qui soient similaires en
# surface. Si vous ne le faites pas, vos blocs aux hautes latitudes seront
# plus petites que les blocs aux basses latitudes.

# Projection Europe
# https://spatialreference.org/ref/esri/europe-lambert-conformal-conic/
projection.europe <- "ESRI:102014"


baseline <- rast("data/baseline.tif")
baseline <- project(baseline,
                    y = projection.europe)
writeRaster(baseline, 
            file = "baseline_lcc.tif",
            overwrite = TRUE)

##### Partie 2 : définition de la taille des blocs ######
# testblocks <-  cv_spatial_autocor(baseline[[c("bio5", "bio6", "forests")]])

# Exploration graphique de la taille des blocs avec shiny :
# cv_block_size(baseline[[ c("bio5", "bio6", "forests")]])

# Identifiez la taille des blocs ici
# Normalement il faudrait prendre le range optimal pour sortir du range 
# d'autocorrélation, mais dans notre cas ça ne marche pas, car ça
# revient à n'avoir qu'une seule cellule ce qui n'est pas envisageable.
# Il faut donc choisir un compromis entre l
#  --> des blocs les plus larges possibles (pour réduire au maximum
# l'autocorrélation entre blocs)
#  --> des blocs suffisamment petits pour pouvoir répartir nos occurrence en
# plusieurs plis (folds)

# Il faut donc procéder par essai-erreur, en créant des blocs larges au début,
# en lançant la boucle de répartition en folds, et en vérifiant le nombre de 
# présences/absences/background par fold. Si c'est équilibré, c'est bon,
# sinon, il faut réessayer avec des blocs plus petits.

# Ici nous allons utiliser le range minimal entre les variables (i.e. range
# d'autocorrélation de la variable forets) :
blocksize <- min(testblocks$range_table$range)
blocksize <- 1000000
# Adaptez cette valeur à vos données !!!!


##### Partie 3: Sélection des backgrounds avant de créer les folds #####

baseline <- rast("data/baseline.tif") 
# Pour les espèces en présence seule : il va falloir générer des folds
# équilibrés en nombre de pseudoabsences ou backgrounds
for (sp in sp_list$sp[which(sp_list$pa.generation == "yes")])
{
  # Points d'occurrence de l'espèce
  P_points <- readRDS(paste0("./data/occurrences_", sp))  
  
  # Vérification pour éviter les accidents bêtes, comme relancer plusieurs fois
  # la boucle...
  if(any(is.na(P_points$Observed)))
  {
    stop("There are NAs in the loaded occurrences. Have you launched the loop",
         "a second time by mistake? If this is intentional, comment/delete the",
         "if() statement that generated this error.")
  }
  
  # Sélectionnez le nombre de backgrounds (plus il y en a, mieux c'est)
  nb_PA <- 1000 
  # Sélectionnez le nombre de répétitions. Notez : si vous avez plus de 10000
  # backgrounds, il n'est probablement pas nécessaire de faire plusieurs 
  # répétitions
  runs_PA <- 3 
  
  
  # Création de l'objet qui contiendra les présences et backgrounds
  P_pa_points <- P_points
  # Cet objet contient actuellement les coordonnées x et y et l'info de présence
  # dans la colonne Observed.
  # Nous allons le compléter en ajoutant, à la suite, les coordonnées des
  # backgrounds pour tous les runs de background
  
  # Pour que biomod sache quelles backgrounds utiliser pour quel run,
  # nous allons créer une deuxième table qui lui dira les lignes correspondant à
  # chaque run de background
  
  # Visualisez cet exercice en écrivant les 2 tableaux en parallèle 
  # (cf. powerpoint)
  
  
  # Création de la matrice d'indication des runs de background 
  # au format attendu par biomod
  bg_biomod_table <- matrix(
    FALSE, # La matrice sera rempli de FALSE par défaut, et on indiquera TRUE 
    # pour les lignes à utiliser
    nrow = nrow(P_points) + runs_PA * nb_PA, # Lignes = présences et backgrounds
    ncol = runs_PA) # Colonnes = runs de background
  # On précise les noms des colonnes au format souhaité par biomod
  colnames(bg_biomod_table) <- paste0("PA", 1:runs_PA)
  
  
  # Quel que soit le run de background, on utilisera toujours les 
  # présences de l'espèce
  # Donc pour tous les runs de background, on met TRUE pour toutes les
  # lignes de présences
  bg_biomod_table[1:nrow(P_points), ] <- TRUE
  
  
  
  # Cette boucle va générer les backgrounds aléatoirement avec la fonction
  # randomPoints du package dismo
  for(PA in 1:runs_PA)
  {
    # D'abord on ajoute les coordonnées des backgrounds à la suite de notre 
    # table de présences
    P_pa_points <- rbind(
      P_pa_points,
      data.frame(spatSample(baseline,
                            size = nb_PA,
                            replace = FALSE, # Pas de remise 
                            na.rm = TRUE, # Pas dans les données manquantes
                            xy = TRUE, # L'output inclut les coords XY
                            values = FALSE), # L'output exclut les variables
                 Observed = NA))
    # Ensuite, on complète notre matrice d'indication des runs de background
    # en mettant TRUE pour toutes les lignes qui correspondent au run qu'on 
    # vient de générer
    bg_biomod_table[(nrow(P_points) + 1 + (PA - 1) * nb_PA):
                             (nrow(P_points) + PA * nb_PA), PA] <- TRUE
  }
  # On sauve les occurrences avec le suffixe "bg" pour dire qu'il y a les
  # points de background dedans
  saveRDS(P_pa_points,
          paste0("data/occurrences_bg_", sp))
  saveRDS(bg_biomod_table,
          paste0("data/bg_table_", sp))
}

##### Etape 4 : préparation des plis (folds) ######
# On va maintenant répartir les blocs en différents plis (=folds), chaque
# pli étant constitué d'un ensemble de blocs qui vont servir à la calibration et
# d'un ensemble de blocs qui vont servir à l'évaluation

# Choisir le nombre de folds (k) de cross-validation
nb_of_folds = 3
# Ratio calibration = (k-1) / k
# Ratio évaluation = 1 / k

# Par exemple, si k = 3 :
# Ratio calibration = 2/3
# Ratio évaluation = 1/3


# Définition des folds
for(i in 1:nrow(sp_list))
{
  sp <- sp_list$sp[i]
  
  # Chargement des occurrences, avec ou sans les backgrounds
  if(sp_list$pa.generation[i] == "yes") {
    P_points <- readRDS(paste0("data/occurrences_bg_", sp))
  } else {
    P_points <- readRDS(paste0("data/occurrences_", sp))
  }
  
  # Transformation des occurrences en objet spatial (format sf)
  pa_data <- sf::st_as_sf(P_points, coords = c("x", "y"), 
                          crs = raster::crs(baseline))
  
  
  if(sp_list$pa.generation[i] == "no") { 
    # Si pas de backgrounds, dans ce cas on ne crée les folds qu'une fois
    cur_folds <- cv_spatial(
      x = pa_data, # Données présence-absence
      column = "Observed", # Colonne contenant l'info présence/absence
      r = baseline[["bio5"]], # Facultatif, quel fond de carte
      # utiliser pour le graphe
      size = blocksize, # Taille des blocs
      k = nb_of_folds, # Nombre de folds
      selection = "random", # Méthode d'attribution des blocs en folds.
      # Cf. ?spatialBlock
      iteration = 50, # Nombre de runs pour répartir les blocs en folds
      biomod2 = TRUE) # Créer une table au format biomod2 ici
    
    # On ne va sauvegarder que la table des folds pour biomod ici
    folds <- cur_folds$biomod_table
    # On renomme les colonnes sinon biomod ne va pas aimer
    colnames(folds) <- paste0("_allData_", 
                              colnames(folds))
    
    # On sauve notre table sur le disque
    saveRDS(folds, paste0("data/folds_", sp, ".RDS"))
    
    png(paste0("./graphiques/folds_", sp, ".png"),
        height = 500, width = 1500)
    p <- cv_plot(cur_folds,
                 x = pa_data,
                 r = baseline[["bio5"]])
    print(p)
    dev.off()
  } else
  {
    # Si backgrounds, dans ce cas on va lancer la procédure de création des folds
    # autant de fois que de runs de PA
    # Chargement de la table des runs backgrounds
    bg_biomod_table <- readRDS(paste0("data/bg_table_", sp))
    
    # Préparation de l'array contenant ligne de calibration et de validation
    # Cf. powerpoint pour mieux comprendre sa conception
    calib_table <- 
      matrix(NA, 
             nrow = nrow(bg_biomod_table),
             ncol = runs_PA * nb_of_folds,
             # Nommer les colonnes à la méthode biomod est un peu technique
             # Il faut avoir le template "_PA1_RUN1" "_PA1_RUN2" "_PA2_RUN1" 
             # etc. On va utiliser sapply() pour le faire correctement
             dimnames = list(NULL,
                             paste0(sapply(paste0("_", 
                                                  colnames(bg_biomod_table)),
                                           rep,
                                           nb_of_folds),
                                    paste0("_RUN", 1:nb_of_folds))
                             )
             )
      

    

    
    for(pa in colnames(bg_biomod_table))
    {
      # On ne sélectionne que les présences et pseuoabsences correspondant au
      # run actuel de backgrounds
      cur_pa_data <- pa_data[which(bg_biomod_table[, pa]), ]
      # Transformation des backgrounds en 0 pour obtenir des folds équilibrés
      # (autant de présences et de backgrounds dans chaque fold)
      cur_pa_data$Observed[which(is.na(cur_pa_data$Observed))] <- 0
      # Si vous souhaitez juste équilibrer les présences, commentez la 
      # ligne au-dessus
      cur_folds <- cv_spatial(
          x = cur_pa_data, # Données présence-absence
          column = "Observed", # Colonne contenant l'info présence/absence
          r = baseline[["forests"]], # Facultatif, quel fond de carte
          # utiliser pour le graphe
          size = blocksize, # Taille des blocs
          k = nb_of_folds, # Nombre de folds
          selection = "random", # Méthode d'attribution des blocs en folds.
          # Cf. ?spatialBlock
          iteration = 50, # Nombre de runs pour répartir les blocs en folds
          biomod2 = TRUE) # Créer une table au format biomod2 ici

      # On remplit notre table
      calib_table[which(bg_biomod_table[, pa]), # Filtrer les lignes à notre run
                  # actuel de PA
                  grep(pa, colnames(calib_table)) # Filtrer les colonnes à notre
                  # run actuel de PA
      ] <- cur_folds$biomod_table
      
      
      png(paste0("graphiques/folds_", sp, "_", pa, ".png"),
          height = 500, width = 1500)
      p <- cv_plot(cur_folds,
                   x = cur_pa_data,
                   r = baseline[["bio5"]])
      print(p)
      dev.off()
    }
    saveRDS(calib_table, 
            paste0("./data/folds_", sp, ".RDS"))
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


library(blockCV)
library(terra)

sp_list <- read.csv("./data_cours/species_list.csv", sep = ";")


##### Partie 1 : Préparation des données #####
# 
# La validation croisée par bloc nécessite d'évaluer le range d'autocorrélation
# des variables environnementales. Pour cela, il faut projeter les raster,
# sinnon le calcul de l'autocorrélation spatiale est imprécis.
# Si vos données ne sont pas projetées (e.g. elles sont en WGS84), il faudra
# les projeter au moins pour le calcul du range d'autocorrélation.
# 
# Si vous travaillez sur de larges étendues spatiales couvrant les hautes 
# latitudes, il vaut mieux projeter toutes vos données à l'étape de 
# préparation des données, afin d'obtenir des blocs qui soient similaires en
# surface. Si vous ne le faites pas, vos blocs aux hautes latitudes seront
# plus petites que les blocs aux basses latitudes.

# Projection Europe
# https://spatialreference.org/ref/esri/europe-lambert-conformal-conic/
projection.europe <- "ESRI:102014"


baseline <- rast("data/baseline.tif")
baseline <- project(baseline,
                    y = projection.europe)
writeRaster(baseline, 
            file = "baseline_lcc.tif",
            overwrite = TRUE)

##### Partie 2 : définition de la taille des blocs ######
# testblocks <-  cv_spatial_autocor(baseline[[c("bio5", "bio6", "forests")]])

# Exploration graphique de la taille des blocs avec shiny :
# cv_block_size(baseline[[ c("bio5", "bio6", "forests")]])

# Identifiez la taille des blocs ici
# Normalement il faudrait prendre le range optimal pour sortir du range 
# d'autocorrélation, mais dans notre cas ça ne marche pas, car ça
# revient à n'avoir qu'une seule cellule ce qui n'est pas envisageable.
# Il faut donc choisir un compromis entre l
#  --> des blocs les plus larges possibles (pour réduire au maximum
# l'autocorrélation entre blocs)
#  --> des blocs suffisamment petits pour pouvoir répartir nos occurrence en
# plusieurs plis (folds)

# Il faut donc procéder par essai-erreur, en créant des blocs larges au début,
# en lançant la boucle de répartition en folds, et en vérifiant le nombre de 
# présences/absences/background par fold. Si c'est équilibré, c'est bon,
# sinon, il faut réessayer avec des blocs plus petits.

# Ici nous allons utiliser le range minimal entre les variables (i.e. range
# d'autocorrélation de la variable forets) :
blocksize <- min(testblocks$range_table$range)
blocksize <- 1000000
# Adaptez cette valeur à vos données !!!!


##### Partie 3: Sélection des backgrounds avant de créer les folds #####

baseline <- rast("data/baseline.tif") 
# Pour les espèces en présence seule : il va falloir générer des folds
# équilibrés en nombre de pseudoabsences ou backgrounds
for (sp in sp_list$sp[which(sp_list$pa.generation == "yes")])
{
  # Points d'occurrence de l'espèce
  P_points <- readRDS(paste0("./data/occurrences_", sp))  
  
  # Vérification pour éviter les accidents bêtes, comme relancer plusieurs fois
  # la boucle...
  if(any(is.na(P_points$Observed)))
  {
    stop("There are NAs in the loaded occurrences. Have you launched the loop",
         "a second time by mistake? If this is intentional, comment/delete the",
         "if() statement that generated this error.")
  }
  
  # Sélectionnez le nombre de backgrounds (plus il y en a, mieux c'est)
  nb_PA <- 1000 
  # Sélectionnez le nombre de répétitions. Notez : si vous avez plus de 10000
  # backgrounds, il n'est probablement pas nécessaire de faire plusieurs 
  # répétitions
  runs_PA <- 3 
  
  
  # Création de l'objet qui contiendra les présences et backgrounds
  P_pa_points <- P_points
  # Cet objet contient actuellement les coordonnées x et y et l'info de présence
  # dans la colonne Observed.
  # Nous allons le compléter en ajoutant, à la suite, les coordonnées des
  # backgrounds pour tous les runs de background
  
  # Pour que biomod sache quelles backgrounds utiliser pour quel run,
  # nous allons créer une deuxième table qui lui dira les lignes correspondant à
  # chaque run de background
  
  # Visualisez cet exercice en écrivant les 2 tableaux en parallèle 
  # (cf. powerpoint)
  
  
  # Création de la matrice d'indication des runs de background 
  # au format attendu par biomod
  bg_biomod_table <- matrix(
    FALSE, # La matrice sera rempli de FALSE par défaut, et on indiquera TRUE 
    # pour les lignes à utiliser
    nrow = nrow(P_points) + runs_PA * nb_PA, # Lignes = présences et backgrounds
    ncol = runs_PA) # Colonnes = runs de background
  # On précise les noms des colonnes au format souhaité par biomod
  colnames(bg_biomod_table) <- paste0("PA", 1:runs_PA)
  
  
  # Quel que soit le run de background, on utilisera toujours les 
  # présences de l'espèce
  # Donc pour tous les runs de background, on met TRUE pour toutes les
  # lignes de présences
  bg_biomod_table[1:nrow(P_points), ] <- TRUE
  
  
  
  # Cette boucle va générer les backgrounds aléatoirement avec la fonction
  # randomPoints du package dismo
  for(PA in 1:runs_PA)
  {
    # D'abord on ajoute les coordonnées des backgrounds à la suite de notre 
    # table de présences
    P_pa_points <- rbind(
      P_pa_points,
      data.frame(spatSample(baseline,
                            size = nb_PA,
                            replace = FALSE, # Pas de remise 
                            na.rm = TRUE, # Pas dans les données manquantes
                            xy = TRUE, # L'output inclut les coords XY
                            values = FALSE), # L'output exclut les variables
                 Observed = NA))
    # Ensuite, on complète notre matrice d'indication des runs de background
    # en mettant TRUE pour toutes les lignes qui correspondent au run qu'on 
    # vient de générer
    bg_biomod_table[(nrow(P_points) + 1 + (PA - 1) * nb_PA):
                             (nrow(P_points) + PA * nb_PA), PA] <- TRUE
  }
  # On sauve les occurrences avec le suffixe "bg" pour dire qu'il y a les
  # points de background dedans
  saveRDS(P_pa_points,
          paste0("data/occurrences_bg_", sp))
  saveRDS(bg_biomod_table,
          paste0("data/bg_table_", sp))
}

##### Etape 4 : préparation des plis (folds) ######
# On va maintenant répartir les blocs en différents plis (=folds), chaque
# pli étant constitué d'un ensemble de blocs qui vont servir à la calibration et
# d'un ensemble de blocs qui vont servir à l'évaluation

# Choisir le nombre de folds (k) de cross-validation
nb_of_folds = 3
# Ratio calibration = (k-1) / k
# Ratio évaluation = 1 / k

# Par exemple, si k = 3 :
# Ratio calibration = 2/3
# Ratio évaluation = 1/3


# Définition des folds
for(i in 1:nrow(sp_list))
{
  sp <- sp_list$sp[i]
  
  # Chargement des occurrences, avec ou sans les backgrounds
  if(sp_list$pa.generation[i] == "yes") {
    P_points <- readRDS(paste0("data/occurrences_bg_", sp))
  } else {
    P_points <- readRDS(paste0("data/occurrences_", sp))
  }
  
  # Transformation des occurrences en objet spatial (format sf)
  pa_data <- sf::st_as_sf(P_points, coords = c("x", "y"), 
                          crs = raster::crs(baseline))
  
  
  if(sp_list$pa.generation[i] == "no") { 
    # Si pas de backgrounds, dans ce cas on ne crée les folds qu'une fois
    cur_folds <- cv_spatial(
      x = pa_data, # Données présence-absence
      column = "Observed", # Colonne contenant l'info présence/absence
      r = baseline[["bio5"]], # Facultatif, quel fond de carte
      # utiliser pour le graphe
      size = blocksize, # Taille des blocs
      k = nb_of_folds, # Nombre de folds
      selection = "random", # Méthode d'attribution des blocs en folds.
      # Cf. ?spatialBlock
      iteration = 50, # Nombre de runs pour répartir les blocs en folds
      biomod2 = TRUE) # Créer une table au format biomod2 ici
    
    # On ne va sauvegarder que la table des folds pour biomod ici
    folds <- cur_folds$biomod_table
    # On renomme les colonnes sinon biomod ne va pas aimer
    colnames(folds) <- paste0("_allData_", 
                              colnames(folds))
    
    # On sauve notre table sur le disque
    saveRDS(folds, paste0("data/folds_", sp, ".RDS"))
    
    png(paste0("./graphiques/folds_", sp, ".png"),
        height = 500, width = 1500)
    p <- cv_plot(cur_folds,
                 x = pa_data,
                 r = baseline[["bio5"]])
    print(p)
    dev.off()
  } else
  {
    # Si backgrounds, dans ce cas on va lancer la procédure de création des folds
    # autant de fois que de runs de PA
    # Chargement de la table des runs backgrounds
    bg_biomod_table <- readRDS(paste0("data/bg_table_", sp))
    
    # Préparation de l'array contenant ligne de calibration et de validation
    # Cf. powerpoint pour mieux comprendre sa conception
    calib_table <- 
      matrix(NA, 
             nrow = nrow(bg_biomod_table),
             ncol = runs_PA * nb_of_folds,
             # Nommer les colonnes à la méthode biomod est un peu technique
             # Il faut avoir le template "_PA1_RUN1" "_PA1_RUN2" "_PA2_RUN1" 
             # etc. On va utiliser sapply() pour le faire correctement
             dimnames = list(NULL,
                             paste0(sapply(paste0("_", 
                                                  colnames(bg_biomod_table)),
                                           rep,
                                           nb_of_folds),
                                    paste0("_RUN", 1:nb_of_folds))
                             )
             )
      

    

    
    for(pa in colnames(bg_biomod_table))
    {
      # On ne sélectionne que les présences et pseuoabsences correspondant au
      # run actuel de backgrounds
      cur_pa_data <- pa_data[which(bg_biomod_table[, pa]), ]
      # Transformation des backgrounds en 0 pour obtenir des folds équilibrés
      # (autant de présences et de backgrounds dans chaque fold)
      cur_pa_data$Observed[which(is.na(cur_pa_data$Observed))] <- 0
      # Si vous souhaitez juste équilibrer les présences, commentez la 
      # ligne au-dessus
      cur_folds <- cv_spatial(
          x = cur_pa_data, # Données présence-absence
          column = "Observed", # Colonne contenant l'info présence/absence
          r = baseline[["forests"]], # Facultatif, quel fond de carte
          # utiliser pour le graphe
          size = blocksize, # Taille des blocs
          k = nb_of_folds, # Nombre de folds
          selection = "random", # Méthode d'attribution des blocs en folds.
          # Cf. ?spatialBlock
          iteration = 50, # Nombre de runs pour répartir les blocs en folds
          biomod2 = TRUE) # Créer une table au format biomod2 ici

      # On remplit notre table
      calib_table[which(bg_biomod_table[, pa]), # Filtrer les lignes à notre run
                  # actuel de PA
                  grep(pa, colnames(calib_table)) # Filtrer les colonnes à notre
                  # run actuel de PA
      ] <- cur_folds$biomod_table
      
      
      png(paste0("graphiques/folds_", sp, "_", pa, ".png"),
          height = 500, width = 1500)
      p <- cv_plot(cur_folds,
                   x = cur_pa_data,
                   r = baseline[["bio5"]])
      print(p)
      dev.off()
    }
    saveRDS(calib_table, 
            paste0("./data/folds_", sp, ".RDS"))
  }
}

>>>>>>> 5768564 (2nd commit)
