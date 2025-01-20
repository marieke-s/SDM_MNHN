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

library(viridis) # Color ramp
library(terra)

#### Section 1 : données environnementales ####
# Pensez à vérifier que les systèmes de coordonnées de vos variables 
# environnementales sont les mêmes, sinon il faut les reprojeter avant
# (fonctions project(), resample())

#---- 1. Données baseline (worldclim 1950-2000) ----
# 1.1 Données climatiques worldclim
climate_vars <- paste0("bio_", 1:19)
env_data <- rast(paste0("data_cours/wc2.1_10m_bio/wc2.1_10m_", climate_vars, 
                        ".tif"))

# Pour avoir les mêmes noms de variables que les fichiers futurs:
names(env_data) <- paste0("bio", 1:19) 

plot(env_data) # Emprise sur le monde entier, par forcément la plus pertinente pour notre espèce

# On réduit l'emprise à une partie de l'ouest-paléarctique
e <- ext(-15, 50, 40, 70)
env_data <- crop(env_data, e)


# Pour tester l'emprise, on peut le dessiner avec la fonction draw :
plot(env_data[[1]])
draw()

# 1.2 Ajout données land use
forests <- rast("data_cours/land use/baseline_landuse.tif")

# Les forêts sont éclatées en plusieurs types de forêts, on les rassemble toutes
# sous une seule variable en faisant la somme 
forests <- sum(forests) 
names(forests) <- "forests"

# Réduction de l'emprise
forests <- crop(forests, e)

# Vérification des résolutions spatiales
res(env_data)
res(forests)

# Les résolutions sont différentes. Stratégie de ré-échantillonnage à définir : 

# Rééchantillonnage : aggregate > resample
# Aggregate possible que si le changement de résolution est un multiple de 2
# Sinon resample
# Comparez
# forests_0.167_agg <- aggregate(x = forests,
#                                fact = 2,
#                                na.rm = TRUE)
# forests_0.167_resample <- resample(x = forests,
#                                    y = env_data)
# plot(c(forests_0.167_resample,
#        forests_0.167_agg))
#
# /!\ Pour le rééchantillonnage, réfléchissez à la méthode à utiliser selon
# les variables
# méthode habituelle pour les variables continues : interpolation bilinéaire
# applicable sur grilles rectilinéaires, interpolation linéaire entre 4 points
# Autres méthodes possibles :
# - Bicubique : version non linéaire, produit une interpolation plus lisse 
#   que la bilinéaire en considérant davantage de points voisins
# - Plus proche voisin : méthode simple, attribue à chaque point la valeur 
#   du point d'échantillonnage le plus proche. Souvent utilisée pour les
#   variables catégorielles (e.g. types de sol, classes de végétation)
# - Inverse Distance Weighting (IDW) : 
#   - Les valeurs sont pondérées par l'inverse de la distance aux points
#     d'échantillonnage voisins
#   - Plus un point est proche, plus il influence la valeur interpolée
#   - Utile pour des variables continues avec une relation spatiale
#     claire (e.g. pollution, température)
# - Interpolation conservative : utilisée pour les variables dont la somme
#   totale doit être conservée après interpolation, e.g. variables 
#   représentant un stock (biomasse, population)
# Si vous travaillez avec des données climatiques de sources variables, 
# il peut être utile de travailler hors de R afin de gérer des formats 
# difficiles à gérer sous R (e.g. données CMIP6), avec python (e.g. outil CDO)

# /!\ Attention : Il ne faut JAMAIS sous-échantillonner car sinon on réplique artificiellement l'information.

forests_0.167 <- resample(forests, 
                          env_data)


env_data <- c(env_data,
              forests_0.167)

# Visualisation
x11() # pour ouvrir une nouvelle fenêtre (uniquement avec Windows)
plot(forests)
x11()
plot(forests_0.167)

# Synchronisation des NAs : explication
val <- values(env_data)
length(which(is.na(val[, "bio1"]) & !is.na(val[, "forests"])))
length(which(is.na(val[, "forests"]) & !is.na(val[, "bio1"])))
# On a des NAs pour différents pixels dans bio et forests 

dummy.raster <- env_data[[1]]
dummy.raster[] <- NA

# Visualisation : 
# Données manquantes pour les variables climatiques = 1
dummy.raster[which(is.na(val[, "bio1"]) & !is.na(val[, "forests"]))] <- 1

# Données manquantes pour les forêts = 2
dummy.raster[which(is.na(val[, "forests"]) & !is.na(val[, "bio1"]))] <- 2

# Nommer les variables catégorielles dans le raster pour faire la carte
levels(dummy.raster) <- data.frame(id = 1:2, label = c("NA climatique", 
                                                       "NA forets"))
plot(dummy.raster, 
     col = c("red", "blue"))

# (Essayez avec plet() pour une version interactive)

# Synchronisation des NAs
# /!\ ATTENTION : il peut être judicieux de faire cette étape plus tard
# Par exemple après la sélection de variables
# Sinon, si vous mettez une variable qui sera plus tard rejetée mais
# qui a beaucoup de NAs, alors elle va rajouter beaucoup de NAs dans les
# autres variables
# Synchronisation 1/3
#
# Fonction mask() : applique un masque (ou gabarit) de valeurs NA sur un autre 
# raster.
# Ici, le masque est défini comme la somme de toutes les variables du raster.
# Si une seule variable contient un NA à un emplacement donné, la somme sera 
# également NA.
# Par conséquent, cet emplacement sera masqué (défini comme NA) dans le raster 
# de sortie.

env_data <- mask(env_data, 
                 app(env_data, fun = sum))

# Ecriture des données
writeRaster(env_data, "data/baseline.tif", overwrite = T) # Overwrite?
rm(list = ls())

#---- 2. Données futures ----
# La difficulté ici est d'associer les variables ensemble par scénario,
# par GCM et par horizon de projection. Les noms diffèrent souvent entre 
# sources, donc il faut regarder
# comment les variables sont nommées sur le disque, et préparer le code pour 
# les lire et les associer dans le bon ordre.
# On va donc commencer par préparer un tableau qui contient les noms des 
# scénarios et GCMs à associer pour chaque source de données (ici : climat et
# land use).

#   -- fichiers climatiques
# GCMs
gcms_clim <- c("IPSL-CM6A-LR", "MIROC6") # GCMs = global climate models

# Scénarios
scenarios_clim<- c("ssp245", "ssp585") # ssp = shared socioeconomic pathways, les scenarios se lisent : ssp 2 4 5, rcp 4.5 ssp = prise en compte des changement géopolitiques, rcp = prise en compte des émissions de gazs à effet de serre.

# Horizons de projection
years_clim <- c("2061-2080")

#   -- fichiers land use
# Notez que l'ordre correspond **exactement** à la section climatique
# GCMs
gcms_lu <- c("ipsl", "miroc") 
# Scénarios
scenarios_lu <- c("rcp45", "rcp85") 
# Horizons de projection
years_lu <- "2070"

# On crée le tableau dans l'ordre scénario > GCM > horizon de projection
proj_names <- data.frame()
for (scenar in 1:length(scenarios_clim)) {
  for (gcm in 1:length(gcms_clim)) {
    for (year in 1:length(years_clim)) {
      proj_names <- rbind.data.frame(proj_names,
                                     data.frame(gcm_clim = gcms_clim[gcm],
                                                gcm_lu = gcms_lu[gcm],
                                                scenar_clim = scenarios_clim[scenar],
                                                scenar_lu = scenarios_lu[scenar],
                                                year_clim = years_clim[year],
                                                year_lu = years_lu[year]))
    }
  }
}

# Vous pouvez aussi préparer ce tableau hors de R, dans excel par exemple,
# si c'est plus simple pour vous.
saveRDS(proj_names, "data/projection_names.RDS")

# On charge les données de baseline
baseline <- rast("data/baseline.tif")

# Important : les noms de variables doivent être identiques entre baseline et 
# futur, on en aura besoin pour projeter les modèles avec biomod

for(i in 1:nrow(proj_names)) {
  # Nom des couches climatiques
  proj_clim <- paste(proj_names[i, c("gcm_clim", "scenar_clim", "year_clim")],
                     collapse = "_")
 
   # Nom des couches land use
  proj_lu <- paste(proj_names[i, c("scenar_lu", "year_lu", "gcm_lu")],
                   collapse = "_")
  
  # Nom finak que l'on souhaite utiliser pour écrire nos couches sur le disque
  # Je choisis un nom simple mais informatif
  final_name <- paste(proj_names[i, c("gcm_lu", "scenar_clim", "year_lu")], 
                      collapse = "_")
  
  cat(final_name, "\n") # imprime le final_name
  
  # Données climatiques
  future_data <- rast(paste0("data_cours/wc2.1_10m_bioc_", 
                             proj_clim,
                             ".tif"))
  names(future_data) <- paste0("bio", 1:19)
  
  # Réduction de l'emprise pour qu'elle corresponde à celle de la baseline
  future_data <- crop(future_data, baseline)
  
  # Données land use
  forests <- rast(paste0("data_cours/land use/LU_", proj_lu, ".tif"))
  forests <- sum(forests) # somme des différents types de forêts
  names(forests) <- "forests" # renommage
  forests <- crop(forests, baseline) # réduction de l'emprise

  forests_0.167 <- resample(forests, 
                            baseline) # rééchantillonnage : utiliser la même méthode partout
  
  future_data <- c(future_data,
                   forests_0.167)
  

  # Synchronisation des NAs entre stacks futurs et synchronisation par rapport 
  # au baseline
  # /!\ ATTENTION : il peut être judicieux de faire cette étape plus tard
  # Par exemple après la sélection de variables
  # Sinon, si vous mettez une variable qui sera plus tard rejetée mais
  # qui a beaucoup de NAs, alors elle va rajouter beaucoup de NAs dans les
  # autres variables
  # Synchronisation 2/3
  
  tmp <- baseline[["bio1"]]
  names(tmp) <- "tmp"
  future_data <- c(tmp, future_data)
  future_data <- mask(future_data, app(future_data, fun = sum))
  future_data <- future_data[[-1]]
  
  writeRaster(future_data, paste("./data/", 
                                 final_name, ".tif", sep = ""), overwrite = T) 
  cat(paste0(Sys.time(), " - ", final_name, "  done\n"))
  
}

# Synchronisation du stack baseline par rapport aux futurs
# /!\ ATTENTION : il peut être judicieux de faire cette étape plus tard
# Par exemple après la sélection de variables
# Sinon, si vous mettez une variable qui sera plus tard rejetée mais
# qui a beaucoup de NAs, alors elle va rajouter beaucoup de NAs dans les
# autres variables
# Synchronisation 3/3
tmp <- future_data[["bio1"]]
names(tmp) <- "tmp"
baseline <- c(tmp,
              baseline)
baseline <- mask(baseline, app(baseline, fun = sum))
baseline <- baseline[[-1]]

writeRaster(baseline, 
            "data/baseline.tif", overwrite = T)

png("graphiques/diff_bio1.png", width = 600, height = 600)
plot(future_data[["bio1"]] - baseline[["bio1"]])
dev.off()



#### Section 2 : données espèce ####
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")

baseline <- rast("data/baseline.tif")
for (sp in sp_list$sp) {
  sp_occurrence <- read.table(paste0("data_cours/", sp, ".csv"), 
                              sep = ";", h = T)
  
  # Pensez à vérifier que les systèmes de coordonnées de vos occurrences
  # correspondent à vos variables environnementales
  # Sinon il faut projeter occurrences ou variables pour assurer la 
  # compatibilité
  
  plot(baseline[[1]])
  points(sp_occurrence[, c("x", "y")], cex = .5, 
         pch = c(1, 16)[sp_occurrence$Observed + 1])
  
  
  # 1. Rasterisation des occurrences
  # D'abord on transforme les occurrences en objet vectoriel spatial
  sp_occurrence <- vect(sp_occurrence,
                        geom = c("x", "y"),
                        crs = "EPSG:4326")
  
  # On transforme le champ "Observed" (composé de 1 & 0) en somme par pixel
  sp_env <- rasterize(x = sp_occurrence, 
                      y = baseline,
                      field = "Observed",
                      fun = sum) # Bien vérifier le résultat et ajuster la 
                                 # fonction
  # Graphe illustrant le nombre de présences par pixel
  plot(sp_env, col = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000"),
       colNA = "lightblue")
  # On ne garde qu'une seule valeur par pixel
  sp_env[sp_env > 1] <- 1
  # On ajoute le nom de l'espèce au raster
  names(sp_env) <- sp
  
  # On ajoute les données environnementales aux présences rasterisées
  sp_env <- c(sp_env,
              baseline)

  # 2. Suppression des présences pour lesquelles on n'a pas de données 
  # environnementales
  # Récupération des coordonnées de toutes les cellules
  coorXY <- xyFromCell(baseline, 1:ncell(baseline))
  # Transformation du raster en data.frame pour tester les NAs
  # Note : cette opération peut être compliquée si vos données environnementales
  # sont trop lourdes, auquel cas il faudra adapter le code
  sp_env_df <- values(sp_env)
  
  # Introduction volontaire d'erreurs pour l'exemple
  # sp_env_df[which(is.na(sp_env_df[, "bio1"])), 1] <- 1
  
  if(any(is.na(sp_env_df[, "bio1"]) & !is.na(sp_env_df[, sp])))
  {
    cat("Some points are in pixels without environmental values\n")
  }
  

  # On supprime les cellules pour lesquelles on n'a pas de données 
  # environnementales
  
  coorXY <- coorXY[-which(is.na(sp_env_df[, "bio1"])), ]
  sp_env_df <- sp_env_df[-which(is.na(sp_env_df[, "bio1"])), ]
  
  # Nombre de cellules résultantes :
  cat(sp, "\nNumber of pixels of presence:",
      "\n - Initial: ", length(which(sp_occurrence$Observed == 1)),
      "\n - After rasterisation: ", length(which(sp_env_df[, 1] == 1)), "\n")
  cat(sp, "\nNumber of pixels of absence:",
      "\n - Initial: ", length(which(sp_occurrence$Observed == 0)),
      "\n - After rasterisation: ", length(which(sp_env_df[, 1] == 0)), "\n\n")
  
  
  # 3. Récupération des occurrences rasterisées et écriture sur le disque
  P_points <- data.frame(
    # D'abord on récupère les coordonnées XY qui correspondent à nos cellules 
    # de présences/absences
    coorXY[which(!is.na(sp_env_df[, sp])), ],
    # Ensuite, on récupère la colonne qui indique présence/absence pour 
    # chaque cellule
    Observed = sp_env_df[which(!is.na(sp_env_df[, sp])), sp]) # On récupère les 
    # occurrences ici

  saveRDS(P_points, file = paste0("data/occurrences_", sp))
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

# Libraries
library(terra)
library(viridis) # Color ramp

#### Section 1 : données environnementales ####
# Pensez à vérifier que les systèmes de coordonnées de vos variables 
# environnementales sont les mêmes, sinon il faut les reprojeter avant
# (fonctions project(), resample())

# 1. Données baseline (worldclim 1950-2000)
# 1.1 Données climatiques worldclim
climate_vars <- paste0("bio_", 1:19)

# Load all climate variables and concatenate into a single raster stack
env_data <- rast(paste0("data_cours/wc2.1_10m_bio/wc2.1_10m_", climate_vars, 
                        ".tif"))

# Pour avoir les mêmes noms de variables que les fichiers futurs:
names(env_data) <- paste0("bio", 1:19) 

plot(env_data) # Etendu spatiale de la variable : le monde entier. Or pour notre espèce ce n'est pas forcément l'étendue la plus pertinente.

# On réduit donc l'emprise à une partie de l'ouest-paléarctique
e <- ext(-15, 50, 40, 70)
env_data <- crop(env_data, e)

# Pour tester l'emprise, on peut le dessiner avec la fonction draw :
plot(env_data[[1]])
draw()


# 1.2 Ajout données land use
forests <- rast("data_cours/land use/baseline_landuse.tif") # baseline = fichier de couverture forestière présent
plot(forests)


# Les forêts sont éclatées en plusieurs types de forêts, on les rassemble toutes sous une seule variable en faisant la somme 
forests <- sum(forests)
plot(forests)

# Réduction de l'emprise
names(forests) <- "forests"
forests <- crop(forests, e)

plot(forests)

# On vérifie les résolutions spatiales des données
res(env_data) #  0.1666667 0.1666667
res(forests) # 0.08333333 0.08333333

# Les résolutions sont différentes 

# Rééchantillonnage : aggregate > resample
# Aggregate possible que si le changement de résolution est un multiple de 2 (ie les grilles sont calées) 
# Sinon resample (si les grilles sont décalées)
# Comparez
# forests_0.167_agg <- aggregate(x = forests,
#                                fact = 2,
#                                na.rm = TRUE)
# forests_0.167_resample <- resample(x = forests,
#                                    y = env_data)
# plot(c(forests_0.167_resample,
#        forests_0.167_agg))
#
# /!\ Pour le rééchantillonnage, réfléchissez à la méthode à utiliser selon
# les variables
# méthode habituelle pour les variables continues : interpolation bilinéaire
# applicable sur grilles rectilinéaires, interpolation linéaire entre 4 points
# Autres méthodes possibles :
# - Bicubique : version non linéaire, produit une interpolation plus lisse 
#   que la bilinéaire en considérant davantage de points voisins
# - Plus proche voisin : méthode simple, attribue à chaque point la valeur 
#   du point d'échantillonnage le plus proche. Souvent utilisée pour les
#   variables catégorielles (e.g. types de sol, classes de végétation)
# - Inverse Distance Weighting (IDW) : 
#   - Les valeurs sont pondérées par l'inverse de la distance aux points
#     d'échantillonnage voisins
#   - Plus un point est proche, plus il influence la valeur interpolée
#   - Utile pour des variables continues avec une relation spatiale
#     claire (e.g. pollution, température)
# - Interpolation conservative : utilisée pour les variables dont la somme
#   totale doit être conservée après interpolation, e.g. variables 
#   représentant un stock (biomasse, population)
# Si vous travaillez avec des données climatiques de sources variables, 
# il peut être utile de travailler hors de R afin de gérer des formats 
# difficiles à gérer sous R (e.g. données CMIP6), avec python (e.g. outil CDO)

# /!\ Attention : Il ne faut JAMAIS sous-échantillonner car sinon on réplique artificiellement l'information 

forests_0.167 <- resample(forests,env_data) # Méthode bilinéaire par défault

# Visualisation
x11() # pour ouvrir une nouvelle fenêtre (uniquement avec Windows)
plot(forests)
x11()
plot(forests_0.167)

# Assemblage des 2 grilles
env_data <- c(env_data,
              forests_0.167)

names(env_data) # We see that we have 20 layers with the last one being "forests"

# Synchronisation des NAs : explication : répercute les NAs présente dans une couche sur toutes les autres couches. 

# Convert raster into df (each column is a layer and each line a grid cell line)
val <- values(env_data)
which(is.na(val[, "bio1"]) & !is.na(val[, "forests"])) # checks NAs in bio1 but not in forests
length(which(is.na(val[, "bio1"]) & !is.na(val[, "forests"]))) # 17 pixels
length(which(is.na(val[, "forests"]) & !is.na(val[, "bio1"]))) # checks number of NAs in forests but not in bio1 : 1699
# bio1 and forests both have different NAs pixels.

dummy.raster <- env_data[[1]]
dummy.raster[] <- NA

# Données manquantes pour les variables climatiques = 1
dummy.raster[which(is.na(val[, "bio1"]) & !is.na(val[, "forests"]))] <- 1

# Données manquantes pour les forêts = 2
dummy.raster[which(is.na(val[, "forests"]) & !is.na(val[, "bio1"]))] <- 2

# Nommer les variables catégorielles dans le raster pour faire la carte
levels(dummy.raster) <- data.frame(id = 1:2, label = c("NA climatique", 
                                                       "NA forets"))
plot(dummy.raster, 
     col = c("red", "blue"))

# (Essayez avec plet() pour une version interactive)

# Synchronisation des NAs
# /!\ ATTENTION : il peut être judicieux de faire cette étape plus tard
# Par exemple après la sélection de variables
# Sinon, si vous mettez une variable qui sera plus tard rejetée mais
# qui a beaucoup de NAs, alors elle va rajouter beaucoup de NAs dans les
# autres variables
# Synchronisation 1/3


# Fonction mask() : applique un masque (ou gabarit) de valeurs NA sur un autre 
# raster.
# Ici, le masque est défini comme la somme de toutes les variables du raster.
# Si une seule variable contient un NA à un emplacement donné, la somme sera 
# également NA.
# Par conséquent, cet emplacement sera masqué (défini comme NA) dans le raster 
# de sortie.
env_data <- mask(env_data, 
                 app(env_data, fun = sum))  # app() applique une fonction à chaque couche du raster


# Ecriture des données
writeRaster(env_data, "data/baseline.tif", overwrite = T) # overwrite = T pour écraser le fichier existant
rm(list = ls())





# 2. Données futures
# La difficulté ici est d'associer les variables ensemble par scénario,
# par GCM et par horizon de projection. Les noms diffèrent souvent entre 
# sources, donc il faut regarder # GCM = global circulation model (Jerome expert de ces modèles)
# comment les variables sont nommées sur le disque, et préparer le code pour 
# les lire et les associer dans le bon ordre.
# On va donc commencer par préparer un tableau qui contient les noms des 
# scénarios et GCMs à associer pour chaque source de données (ici : climat et
# land use).

#   -- fichiers climatiques
# GCMs
gcms_clim <- c("IPSL-CM6A-LR", "MIROC6") # Toujours prendre plusieurs modèles

# Scénarios
scenarios_clim<- c("ssp245", "ssp585") # scénarios climatiques (ssp : socio-economic pathways - scenarios du dernier rapport des experts climatiques - prend en compte les changements géopolitiques globaux), se lit ssp 2 rcp 4.5 (scenario d'emission de gazs à effet de serre).

# Horizons de projection
years_clim <- c("2061-2080")

#   -- fichiers land use
# Notez que l'ordre correspond **exactement** à la section climatique
# GCMs
gcms_lu <- c("ipsl", "miroc") 

# Scénarios
scenarios_lu <- c("rcp45", "rcp85") 

# Horizons de projection
years_lu <- "2070"

# On crée le tableau dans l'ordre scénario > GCM > horizon de projection
proj_names <- data.frame()
for (scenar in 1:length(scenarios_clim)) {
  for (gcm in 1:length(gcms_clim)) {
    for (year in 1:length(years_clim)) {
      proj_names <- rbind.data.frame(proj_names,
                                     data.frame(gcm_clim = gcms_clim[gcm],
                                                gcm_lu = gcms_lu[gcm],
                                                scenar_clim = scenarios_clim[scenar],
                                                scenar_lu = scenarios_lu[scenar],
                                                year_clim = years_clim[year],
                                                year_lu = years_lu[year]))
    }
  }
}

# Vous pouvez aussi préparer ce tableau hors de R, dans excel par exemple, si c'est plus simple pour vous.
saveRDS(proj_names, "data/projection_names.RDS") # RDS 100 x plus compressé que CSV, ne fonctionne pas pour les raster terra (qui ne sont pas chargés dans la mémoire)

# On charge les données de base
baseline <- rast("data/baseline.tif")

# Important : les noms de variables doivent être identiques entre baseline et 
# futur, on en aura besoin pour projeter les modèles avec biomod
for(i in 1:nrow(proj_names)) {
  # Nom des couches climatiques
  proj_clim <- paste(proj_names[i, c("gcm_clim", "scenar_clim", "year_clim")],
                     collapse = "_")
  # Nom des couches land use
  proj_lu <- paste(proj_names[i, c("scenar_lu", "year_lu", "gcm_lu")],
                   collapse = "_")
  
  # Nom finak que l'on souhaite utiliser pour écrire nos couches sur le disque
  # Je choisis un nom simple mais informatif
  final_name <- paste(proj_names[i, c("gcm_lu", "scenar_clim", "year_lu")], 
                      collapse = "_")
  
  cat(final_name, "\n")
  
  # Données climatiques
  future_data <- rast(paste0("data_cours/wc2.1_10m_bioc_", 
                             proj_clim,
                             ".tif"))
  names(future_data) <- paste0("bio", 1:19)
  future_data <- crop(future_data, baseline)
  
  
  # Données land use
  forests <- rast(paste0("data_cours/land use/LU_", proj_lu, ".tif"))
  forests <- sum(forests)
  names(forests) <- "forests"
  forests <- crop(forests, baseline)

  forests_0.167 <- resample(forests, 
                            baseline)
  
  future_data <- c(future_data,
                   forests_0.167)
  

  # Synchronisation des NAs entre stacks futurs et synchronisation par rapport 
  # au baseline
  # /!\ ATTENTION : il peut être judicieux de faire cette étape plus tard
  # Par exemple après la sélection de variables
  # Sinon, si vous mettez une variable qui sera plus tard rejetée mais
  # qui a beaucoup de NAs, alors elle va rajouter beaucoup de NAs dans les
  # autres variables
  # Synchronisation 2/3
  
  tmp <- baseline[["bio1"]]
  names(tmp) <- "tmp"
  future_data <- c(tmp, future_data)
  future_data <- mask(future_data, app(future_data, fun = sum))
  future_data <- future_data[[-1]]
  
  writeRaster(future_data, paste("./data/", 
                                 final_name, ".tif", sep = ""), overwrite = T) 
  cat(paste0(Sys.time(), " - ", final_name, "  done\n"))
  
}

# Synchronisation du stack baseline par rapport aux futurs
# /!\ ATTENTION : il peut être judicieux de faire cette étape plus tard
# Par exemple après la sélection de variables
# Sinon, si vous mettez une variable qui sera plus tard rejetée mais
# qui a beaucoup de NAs, alors elle va rajouter beaucoup de NAs dans les
# autres variables
# Synchronisation 3/3
tmp <- future_data[["bio1"]]
names(tmp) <- "tmp"
baseline <- c(tmp,
              baseline)
baseline <- mask(baseline, app(baseline, fun = sum))
baseline <- baseline[[-1]]

writeRaster(baseline, 
            "data/baseline.tif", overwrite = T)

png("graphiques/diff_bio1.png", width = 600, height = 600)
plot(future_data[["bio1"]] - baseline[["bio1"]])
dev.off()



#### Section 2 : données espèce ####
sp_list <- read.csv("data_cours/species_list.csv", sep = ";")

baseline <- rast("data/baseline.tif")
for (sp in sp_list$sp) {
  sp_occurrence <- read.table(paste0("data_cours/", sp, ".csv"), 
                              sep = ";", h = T)
  
  # Pensez à vérifier que les systèmes de coordonnées de vos occurrences correspondent à vos variables environnementales
  # Sinon il faut projeter occurrences ou variables pour assurer la compatibilité
  
  plot(baseline[[1]])
  points(sp_occurrence[, c("x", "y")], cex = .5, 
         pch = c(1, 16)[sp_occurrence$Observed + 1])
  
  
  # 1. Rasterisation des occurrences
  # D'abord on transforme les occurrences en objet vectoriel spatial
  sp_occurrence <- vect(sp_occurrence,
                        geom = c("x", "y"),
                        crs = "EPSG:4326")
  
  # On transforme le champ "Observed" (composé de 1 & 0) en somme par pixel
  sp_env <- rasterize(x = sp_occurrence, 
                      y = baseline,
                      field = "Observed",
                      fun = sum) # Bien vérifier le résultat et ajuster la 
                                 # fonction
  # Graphe illustrant le nombre de présences par pixel
  plot(sp_env, col = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000"),
       colNA = "lightblue")
  # On ne garde qu'une seule valeur par pixel
  sp_env[sp_env > 1] <- 1
  # On ajoute le nom de l'espèce au raster
  names(sp_env) <- sp
  
  # On ajoute les données environnementales aux présences rasterisées
  sp_env <- c(sp_env,
              baseline)

  # 2. Suppression des présences pour lesquelles on n'a pas de données 
  # environnementales
  # Récupération des coordonnées de toutes les cellules
  coorXY <- xyFromCell(baseline, 1:ncell(baseline))
  # Transformation du raster en data.frame pour tester les NAs
  # Note : cette opération peut être compliquée si vos données environnementales
  # sont trop lourdes, auquel cas il faudra adapter le code
  sp_env_df <- values(sp_env)
  
  # Introduction volontaire d'erreurs pour l'exemple
  # sp_env_df[which(is.na(sp_env_df[, "bio1"])), 1] <- 1
  
  if(any(is.na(sp_env_df[, "bio1"]) & !is.na(sp_env_df[, sp])))
  {
    cat("Some points are in pixels without environmental values\n")
  }
  

  # On supprime les cellules pour lesquelles on n'a pas de données 
  # environnementales
  
  coorXY <- coorXY[-which(is.na(sp_env_df[, "bio1"])), ]
  sp_env_df <- sp_env_df[-which(is.na(sp_env_df[, "bio1"])), ]
  
  # Nombre de cellules résultantes :
  cat(sp, "\nNumber of pixels of presence:",
      "\n - Initial: ", length(which(sp_occurrence$Observed == 1)),
      "\n - After rasterisation: ", length(which(sp_env_df[, 1] == 1)), "\n")
  cat(sp, "\nNumber of pixels of absence:",
      "\n - Initial: ", length(which(sp_occurrence$Observed == 0)),
      "\n - After rasterisation: ", length(which(sp_env_df[, 1] == 0)), "\n\n")
  
  
  # 3. Récupération des occurrences rasterisées et écriture sur le disque
  P_points <- data.frame(
    # D'abord on récupère les coordonnées XY qui correspondent à nos cellules 
    # de présences/absences
    coorXY[which(!is.na(sp_env_df[, sp])), ],
    # Ensuite, on récupère la colonne qui indique présence/absence pour 
    # chaque cellule
    Observed = sp_env_df[which(!is.na(sp_env_df[, sp])), sp]) # On récupère les 
    # occurrences ici

  saveRDS(P_points, file = paste0("data/occurrences_", sp))
}

>>>>>>> 5768564 (2nd commit)
