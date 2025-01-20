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

sp_list <- read.csv("./data_cours/species_list.csv", sep = ";")

proj_names <- readRDS("data/projection_names.RDS")

list_stacks <- c("baseline",
                 apply(proj_names[, c("gcm_lu", "scenar_clim", "year_lu")], 1, 
                       paste, collapse = "_"))

scenarios <- c("ssp245_2070", "ssp585_2070")

for (i in 1:nrow(sp_list))
{
  sp <- sp_list$sp[i]
  cat(paste("----", Sys.time(), sp, "stack creation initialised ----\n", 
            sep = " "))
 
  
  # (Facultatif) On va charger les coordonnées de présence pour faire des cartes 
  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))
  input_data <- get_formal_data(model_runs)
  sp_coords <- input_data@coord[which(input_data@data.species == 1), ] 
  
  # On commence par charger toutes les projections que l'on a faites
  cat(paste("-- Stacking continuous maps...\n", sep = " "))
  for(s in list_stacks)
  {
    cur_em <- rast(paste0("models/", sp, "/proj_", s,
                          "/proj_", s, "_", sp, "_ensemble.tif"))
    if(s == list_stacks[1])
    {
      # Choisir le bon modèle d'ensemble names(cur_em)
      continuous_stack_cm <- cur_em[[grep("EMmean", names(cur_em))]]
    } else
    {
      continuous_stack_cm <- c(continuous_stack_cm,
                               cur_em[[grep("EMmean", names(cur_em))]])
    }
  }
  names(continuous_stack_cm) <- list_stacks

  # continuous_stack_cm contient toutes les projections, et donc pour chaque
  # scénario on a plusieurs projections car il y a plusieurs GCMs.
  # Ce qui nous intéresse c'est d'avoir une seule projection par GCM
  # Donc la prochaine étape va consister à créer un stack qui ne contient
  # qu'une carte par scénario :
  # La première couche est le baseline, cette couche ne change pas
  continuous_stack <- continuous_stack_cm[[1]]
  
  cat(paste("-- Plotting continuous stacks...\n"))
  # En bonus on fait les cartes au format pdf
  pdf(paste0("graphiques/", sp, "_cont.pdf", sep = ""), width = 11.7, 
      height = 8.3)
  par(mar = c(2.1, 2.6, 2.6, 3.6))
  # Carte du baseline
  plot(continuous_stack[[1]],
       las = 1, 
       range = c(0, 1000), 
       main = paste(sp, " baseline", sep = ""),
       col = rev(viridis::inferno(100))) # Requires package viridis
  points(sp_coords, cex = .4, pch = 16, col = "cyan")
  
  # Ensuite on va faire une seule carte par scénario avec une boucle
  # Note : scénario sous-entend ici scénario + horizon de projection temporel
  # e.g. SSP5-8.5 2070
  for (k in scenarios)
  {
    # On ne conserve que les projections (GCMs) qui correspondent à notre 
    # scénario avec grep() qui recherche les mots qu'on lui donne
    subStack <- continuous_stack_cm[[grep(k, names(continuous_stack_cm))]]
    # Pour n'avoir qu'une seule carte par scénario, on fait la moyenne des
    # GCMs du scénario :
    continuous_stack <- c(continuous_stack,
                          mean(subStack))
    
    # On peut alors projeter la carte du scénario
    plot(continuous_stack[[nlyr(continuous_stack)]],
         las = 1, 
         range = c(0, 1000), 
         main = paste(sp, " ", k, sep = ""),
         col = rev(viridis::inferno(100)))
    
    points(sp_coords, cex = .4, pch = 16, col = "cyan")
    
  }
  dev.off()
  names(continuous_stack) <- c("baseline", scenarios)
  cat(paste("-- Writing stacks...\n", sep = " "))
  if(!dir.exists("outputs"))
  {
    dir.create("outputs", recursive = T)
  }
  writeRaster(continuous_stack, paste0("outputs/cont_", sp, ".tif"), 
              overwrite = T)
  cat(paste("----", Sys.time(), sp, "continuous stack creation finished ----\n", 
            sep = " "))
  
  
  cat(paste("-- Stacking binary maps...\n", sep = " "))
  # La suite est quasiment la même opération que pour les stacks continus
  # Sauf que l'on ne fait pas la moyenne des GCMs, on va plutôt faire du
  # committee averaging
  for(s in list_stacks)
  {
    cur_em <- rast(paste0("models/", sp, "/proj_", s,
                          "/proj_", s, "_", sp, "_ensemble_TSSbin.tif"))
    if(s == list_stacks[1])
    {
      binary_stack_cm <- cur_em[[grep("EMmean", names(cur_em))]] 
    } else
    {
      binary_stack_cm <- c(binary_stack_cm,
                           cur_em[[grep("EMmean", names(cur_em))]])
    }
  }
  names(binary_stack_cm) <- list_stacks
  
  binary_stack <- binary_stack_cm[[1]]
  
  cat(paste("-- Plotting binary stacks...\n"))
  pdf(paste0("graphiques/", sp, "_binary.pdf", sep = ""), width = 11.7, height = 8.3)
  par(mar = c(2.1, 2.6, 2.6, 3.6))
  plot(binary_stack[[1]],
       las = 1,  
       main = paste(sp, " baseline", sep = ""))
  points(sp_coords, cex = .5, pch = 16)
  for (k in scenarios)
  {
    subStack <- binary_stack_cm[[grep(k, names(binary_stack_cm))]]
    # Committee averaging (démocratie)
    # Si plus de 50% des projections prédisent présence on met présence, sinon
    # on met absence
    bin_tmp <- sum(subStack)
    bin_tmp[bin_tmp <= nlyr(subStack)/2] <- 0
    bin_tmp[bin_tmp > nlyr(subStack)/2] <- 1
    
    binary_stack <- c(binary_stack,
                      bin_tmp)
    plot(binary_stack[[nlyr(binary_stack)]],
         las = 1,
         main = paste(sp, " ", k, sep = ""))
    
    points(sp_coords, cex = .5, pch = 16)
    
  }
  dev.off()
  names(binary_stack) <- c("baseline", scenarios)
  cat(paste("-- Writing stacks...\n", sep = " "))
  
  writeRaster(binary_stack, paste0("outputs/bin_", sp, ".tif"), overwrite = T)
  cat(paste("----", Sys.time(), sp, "binary stack creation finished ----\n", sep = " "))
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


library(biomod2)

sp_list <- read.csv("./data_cours/species_list.csv", sep = ";")

proj_names <- readRDS("data/projection_names.RDS")

list_stacks <- c("baseline",
                 apply(proj_names[, c("gcm_lu", "scenar_clim", "year_lu")], 1, 
                       paste, collapse = "_"))

scenarios <- c("ssp245_2070", "ssp585_2070")

for (i in 1:nrow(sp_list))
{
  sp <- sp_list$sp[i]
  cat(paste("----", Sys.time(), sp, "stack creation initialised ----\n", 
            sep = " "))
 
  
  # (Facultatif) On va charger les coordonnées de présence pour faire des cartes 
  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))
  input_data <- get_formal_data(model_runs)
  sp_coords <- input_data@coord[which(input_data@data.species == 1), ] 
  
  # On commence par charger toutes les projections que l'on a faites
  cat(paste("-- Stacking continuous maps...\n", sep = " "))
  for(s in list_stacks)
  {
    cur_em <- rast(paste0("models/", sp, "/proj_", s,
                          "/proj_", s, "_", sp, "_ensemble.tif"))
    if(s == list_stacks[1])
    {
      # Choisir le bon modèle d'ensemble names(cur_em)
      continuous_stack_cm <- cur_em[[grep("EMmean", names(cur_em))]]
    } else
    {
      continuous_stack_cm <- c(continuous_stack_cm,
                               cur_em[[grep("EMmean", names(cur_em))]])
    }
  }
  names(continuous_stack_cm) <- list_stacks

  # continuous_stack_cm contient toutes les projections, et donc pour chaque
  # scénario on a plusieurs projections car il y a plusieurs GCMs.
  # Ce qui nous intéresse c'est d'avoir une seule projection par GCM
  # Donc la prochaine étape va consister à créer un stack qui ne contient
  # qu'une carte par scénario :
  # La première couche est le baseline, cette couche ne change pas
  continuous_stack <- continuous_stack_cm[[1]]
  
  cat(paste("-- Plotting continuous stacks...\n"))
  # En bonus on fait les cartes au format pdf
  pdf(paste0("graphiques/", sp, "_cont.pdf", sep = ""), width = 11.7, 
      height = 8.3)
  par(mar = c(2.1, 2.6, 2.6, 3.6))
  # Carte du baseline
  plot(continuous_stack[[1]],
       las = 1, 
       range = c(0, 1000), 
       main = paste(sp, " baseline", sep = ""),
       col = rev(viridis::inferno(100))) # Requires package viridis
  points(sp_coords, cex = .4, pch = 16, col = "cyan")
  
  # Ensuite on va faire une seule carte par scénario avec une boucle
  # Note : scénario sous-entend ici scénario + horizon de projection temporel
  # e.g. SSP5-8.5 2070
  for (k in scenarios)
  {
    # On ne conserve que les projections (GCMs) qui correspondent à notre 
    # scénario avec grep() qui recherche les mots qu'on lui donne
    subStack <- continuous_stack_cm[[grep(k, names(continuous_stack_cm))]]
    # Pour n'avoir qu'une seule carte par scénario, on fait la moyenne des
    # GCMs du scénario :
    continuous_stack <- c(continuous_stack,
                          mean(subStack))
    
    # On peut alors projeter la carte du scénario
    plot(continuous_stack[[nlyr(continuous_stack)]],
         las = 1, 
         range = c(0, 1000), 
         main = paste(sp, " ", k, sep = ""),
         col = rev(viridis::inferno(100)))
    
    points(sp_coords, cex = .4, pch = 16, col = "cyan")
    
  }
  dev.off()
  names(continuous_stack) <- c("baseline", scenarios)
  cat(paste("-- Writing stacks...\n", sep = " "))
  if(!dir.exists("outputs"))
  {
    dir.create("outputs", recursive = T)
  }
  writeRaster(continuous_stack, paste0("outputs/cont_", sp, ".tif"), 
              overwrite = T)
  cat(paste("----", Sys.time(), sp, "continuous stack creation finished ----\n", 
            sep = " "))
  
  
  cat(paste("-- Stacking binary maps...\n", sep = " "))
  # La suite est quasiment la même opération que pour les stacks continus
  # Sauf que l'on ne fait pas la moyenne des GCMs, on va plutôt faire du
  # committee averaging
  for(s in list_stacks)
  {
    cur_em <- rast(paste0("models/", sp, "/proj_", s,
                          "/proj_", s, "_", sp, "_ensemble_TSSbin.tif"))
    if(s == list_stacks[1])
    {
      binary_stack_cm <- cur_em[[grep("EMmean", names(cur_em))]] 
    } else
    {
      binary_stack_cm <- c(binary_stack_cm,
                           cur_em[[grep("EMmean", names(cur_em))]])
    }
  }
  names(binary_stack_cm) <- list_stacks
  
  binary_stack <- binary_stack_cm[[1]]
  
  cat(paste("-- Plotting binary stacks...\n"))
  pdf(paste0("graphiques/", sp, "_binary.pdf", sep = ""), width = 11.7, height = 8.3)
  par(mar = c(2.1, 2.6, 2.6, 3.6))
  plot(binary_stack[[1]],
       las = 1,  
       main = paste(sp, " baseline", sep = ""))
  points(sp_coords, cex = .5, pch = 16)
  for (k in scenarios)
  {
    subStack <- binary_stack_cm[[grep(k, names(binary_stack_cm))]]
    # Committee averaging (démocratie)
    # Si plus de 50% des projections prédisent présence on met présence, sinon
    # on met absence
    bin_tmp <- sum(subStack)
    bin_tmp[bin_tmp <= nlyr(subStack)/2] <- 0
    bin_tmp[bin_tmp > nlyr(subStack)/2] <- 1
    
    binary_stack <- c(binary_stack,
                      bin_tmp)
    plot(binary_stack[[nlyr(binary_stack)]],
         las = 1,
         main = paste(sp, " ", k, sep = ""))
    
    points(sp_coords, cex = .5, pch = 16)
    
  }
  dev.off()
  names(binary_stack) <- c("baseline", scenarios)
  cat(paste("-- Writing stacks...\n", sep = " "))
  
  writeRaster(binary_stack, paste0("outputs/bin_", sp, ".tif"), overwrite = T)
  cat(paste("----", Sys.time(), sp, "binary stack creation finished ----\n", sep = " "))
}
>>>>>>> 5768564 (2nd commit)
