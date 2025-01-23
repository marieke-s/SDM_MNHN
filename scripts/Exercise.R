# Exercice : 
# Calibrer 2 modèles sans biomod et sortir les cartes (pour une seule espèces) :
# 1. RF non paramétré
# 2. RF bien paramétré


#---- Setting up ----
# Load libraries
library(biomod2) # Testé fonctionnel en v. 4.2-5-2 --> pb de rétrocompatibilité avec les versions antérieures possible
library(ggtext)
library(RColorBrewer) # Pour les palettes de couleurs
library(terra)
library(tidyterra)
library(randomForest)

#---- Data prep ----
# Load necessary data 
sp <- "Jon.snowii"
P_points <- readRDS(paste0("./data/occurrences_", sp))
baseline <- rast("data/baseline.tif")

# Ajout des coordonnées des backgrounds à la suite de notre table de présences
P_pa_points <- rbind(
  P_points,
  data.frame(spatSample(baseline,
                        size = 15000,
                        replace = FALSE, # Pas de remise 
                        na.rm = TRUE, # Pas dans les données manquantes
                        xy = TRUE, # L'output inclut les coords XY
                        values = FALSE), # L'output exclut les variables # Si df, on met TRUE ici
             Observed = NA))

# Ajout des données environnementales bio5, bio6 et forests au tableau :
# Extract data 
pred <- terra::extract(baseline[[c("bio5", "bio6", "forests")]], P_pa_points[, c("x", "y")], ID = FALSE)

# Bind PA and pred
P_pa_points <- cbind(P_pa_points, pred)
head(P_pa_points)

# Remplacement des NAs en 0
P_pa_points$Observed[which(is.na(P_pa_points$Observed))] <- 0

# Variable selection 
# Already done here : bio 5, bio 6 and forests

# Block cross-validation
# On ne fait pas ça dans l'exercise

#---- Calibration and Predictionavec le package RF ----

P_pa_points$Observed <- as.factor(P_pa_points$Observed)
rf0 <- randomForest(Observed ~ bio5 + bio6 + forests, data = P_pa_points, norm.votes = TRUE)


# Use terra package to make the prediction

predict(baseline[[c("bio5", "bio6", "forests")]], rf0)

# Error : we need to remove NAs

p <- predict(baseline[[c("bio5", "bio6", "forests")]], rf0, na.rm = TRUE)
plot(p)

# It gives only presence and absence but not the proba : we need to add the parameter type = "prob" from the randomForest::predict arguments.

p <- predict(baseline[[c("bio5", "bio6", "forests")]], rf0, na.rm = TRUE, type = "prob")
plot(p)

# On obtient 2 cartes : X0 la proba d'avoir des absences et X1 la proba d'avoir des présences (X0 = 1- X1)
# Le modèle est overfitté 

# Paramétrisation du RF
# Downsampling :Chaque arbre est calculé sur un nombre équilibré de présence et d'abscence --> 

# Nb de présences: 
prNum <- length(which(P_pa_points$Observed == 1))

# Creation d'un vecteur équilibré: 
samplesize <- c("0" = prNum, "1" = prNum)

# On recalibre le modèle avec le nouveau vecteur
rf1 <- randomForest(Observed ~ bio5 + bio6 + forests, data = P_pa_points, norm.votes = TRUE, sampsize = samplesize)

# On refait la prédiction
p <- predict(baseline[[c("bio5", "bio6", "forests")]], rf1, na.rm = TRUE, type = "prob")
plot(p)

# On obtient des probabilités plus réalistes


