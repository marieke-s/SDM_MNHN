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
library(terra)
library(reshape2)

sp_list <- read.csv("./data_cours/species_list.csv", sep = ";")

for (i in 1:nrow(sp_list)) 
{
  sp <- sp_list$sp[i]
  
  # Chargement des données formatées et des modèles calibrés
  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))
  input_data <- get_formal_data(model_runs)
  
  # Variables utilisées pour la calibration
  cur_vars <- model_runs@expl.var.names
  
  # Calcul des courbes de réponse
  resp <- bm_PlotResponseCurves(bm.out = model_runs,
                                fixed.var = "mean",
                                data_species = input_data@data.species
                                )$tab
  # Faites la comparaison en enlevant le data_species !!
    
  colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")

  # Ordonnons les variables dans l'ordre d'entrée
  resp$Variable <- factor(resp$Variable, levels = cur_vars)
  
  p <- ggplot(resp, aes(x = Var.value, y = Response))+ 
    geom_line(alpha = 0.2, aes(group = Model)) + 
    stat_smooth() +
    facet_wrap(~Variable, scales = "free_x") + 
    theme_bw() + 
    ylim(0, 1.1) + 
    xlab("Variable value")
  
  # A réfléchir : ajouter limites des données de présence ?
  
  png(paste0("graphiques/response_plot_", sp, ".png"), width = 550 * 4.2, height = 550 * 4.2, res = 300)
  
  print(p)
  
  dev.off()
}
