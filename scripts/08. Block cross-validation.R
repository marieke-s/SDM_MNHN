# Clean env and plots 
rm(list = ls())
dev.off()

# Load libraries
library(blockCV)
library(terra)

# Load species list
sp_list <- read.csv("./data_cours/species_list.csv", sep = ";")

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

# Load baseline predictors
baseline <- rast("data/baseline.tif")

# Check de la taille des pixels projetés 
plot(terra::cellSize(baseline))

# Projection Europe
# https://spatialreference.org/ref/esri/europe-lambert-conformal-conic/
projection.europe <- "ESRI:102014"

# Projection du raster 
baseline <- project(baseline,
                    y = projection.europe)

# Check de la taille des pixels projetés 
plot(terra::cellSize(baseline))

# Sauvegarde du raster projeté
writeRaster(baseline, 
            file = "baseline_lcc.tif",
            overwrite = TRUE)


# Attention : quand on projete on change les valeurs des pixels.




#----------- Parenthèse :  Comparaison de la surface d'occurence ----------
# Load raster with predicted occurence 
a <- rast("./outputs/em_bin_Jon.snowii.tif")
plot(a)

plot(a * cellSize(a))

global(a * cellSize(a), sum, na.rm = TRUE)

rm(a)
##### Partie 2 : définition de la taille des blocs ######

# Calcul du range d'autocorrelation (sur nos varaibles finales en entier (mettre les raster complets c'est mieux car ça donne plus de pixels pour faire le calcul plus que de donner que les points d'occurrences.))
testblocks <-  cv_spatial_autocor(baseline[[c("bio5", "bio6", "forests")]])

# On obtient une taille de bloc qui fait plusieurs fois la Terre. 
# Il va falloir tester un seuil plus bas.

# Exploration graphique de la taille des blocs avec shiny :
 cv_block_size(baseline[[ c("bio5", "bio6", "forests")]])
 
 # Chosir la taille la plus grande possible tout en ayant la possibilité de faire plusieurs folds équilibrés.--> c'est sur d'avoir des blocs équilibrés quand ils sont trop grand. 

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
blocksize <- 2000000

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
  # spatSample du package terra
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
                            values = FALSE), # L'output exclut les variables # Si df, on met TRUE ici
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
# le nb de fold determine le nb de points qu'on a dans chaque calib/test set
# Chaque fold produit 1 modèle avec un jeu de calibration (un groupe de blocs) et un jeu de test (groupe de bloc) différents. --> On doit regarder ensuite si les résultats entre les folds sont stables. Inquiétez vous si un des fold a une résultat très différent --> la variabilité peut indiquer que l'exclusion de certains blocs ne permettent pas de bien calibrer votre md. --> notamment si on des conditions assez hétérogènes : par ex si des points sont au nord avec des conditions particulières. 
# Choix du nb de fold : tatonnement, on essaie entre 3 et 5 et puis on regarde si les folds sont équilibrés. Minimum 3 folds, jusqu'à 10 folds.
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
  # Pcq c'est ce que blockCV demande
  pa_data <- sf::st_as_sf(P_points, coords = c("x", "y"), 
                          crs = raster::crs(baseline)) # On retourne en WGS84 mais on aurait pu aussi tout faire en projeté

  
  
  if(sp_list$pa.generation[i] == "no") { 
     # Si pas de backgrounds, dans ce cas on ne crée les folds qu'une fois :
    
    cur_folds <- cv_spatial(
      x = pa_data, # Données présence-absence (objet créé dans ce script au dessus)
      column = "Observed", # Colonne contenant l'info présence/absence
      r = baseline[["bio5"]], # Facultatif, quel fond de carte utiliser pour le graphe
      size = blocksize, # Taille des blocs
      k = nb_of_folds, # Nombre de folds
      selection = "random", # Méthode d'attribution des blocs en folds.Cf. ?spatialBlock
      iteration = 50, # Nombre de runs pour répartir les blocs en folds. Nombre d'essais pour trouver des plis équilibrés --> prend le meilleur, le + équilibré
      biomod2 = TRUE) # Créer une table au format biomod2 ici
    
    # Vérifier le tableau pour voir s'il est déséquilibrer. Si vous avez beaucoup d'espèces il faut implémenter un code qui test automatiquement. Par exemple : 
        # Test if folds are balanced
    # ifelse(max((folds$records$train_1 - median(folds$records$train_1)) / median(folds$records$train_1)) > 0.3, "déséquilibrés", "équilibrés")
    
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
    # autant de fois que de runs de PA (cf tableau ppt):
      
    # Chargement de la table des runs backgrounds
    bg_biomod_table <- readRDS(paste0("data/bg_table_", sp))
    
    # Préparation de l'array contenant ligne de calibration et de validation
    # Cf. powerpoint pour mieux comprendre sa conception
    calib_table <- 
      matrix(NA, 
             nrow = nrow(bg_biomod_table), # nombre de lignes
             ncol = runs_PA * nb_of_folds, # nombre de colonnes
             # Nommer les colonnes à la méthode biomod est un peu technique
             # Il faut avoir le template "_PA1_RUN1" "_PA1_RUN2" "_PA2_RUN1" 
             # etc. On va utiliser sapply() pour le faire correctement
             dimnames = list(NULL, # pas de nom pour les lignes
                             paste0(sapply(paste0("_", 
                                                  colnames(bg_biomod_table)),
                                           rep,
                                           nb_of_folds), # je repète 'PA' autant de fois que j'ai de folds
                                    paste0("_RUN", 1:nb_of_folds))
                             )
             )
      
    
    for(pa in colnames(bg_biomod_table))
   
       {
      
      # On ne sélectionne que les présences et pseuoabsences correspondant au
      # run actuel de backgrounds
      cur_pa_data <- pa_data[which(bg_biomod_table[, pa]), ]
      
      # Transformation des backgrounds (NA) en 0 pour obtenir des folds équilibrés
      # temporairement pour que ça fonctionne avec blockCV
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
          selection = "random", # Méthode d'attribution des blocs en folds. Cf. ?spatialBlock
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


