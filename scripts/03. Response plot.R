
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
library(ggplot2)
library(terra)
library(reshape2)

sp_list <- read.csv("./data_cours/species_list.csv", sep = ";")

for (i in 1:nrow(sp_list)) 
{
  sp <- sp_list$sp[i]
  
  # Chargement des données formatées et des modèles calibrés pour connaitre les coordonnées des spp. pour pouvoir faire la moyenne des variables en ce point.
  model_runs <- readRDS(paste0("models/", sp, "/model_runs.RDS"))
  input_data <- get_formal_data(model_runs) # récupération des données d'entrée du md.
  
  # Variables utilisées pour la calibration
  cur_vars <- model_runs@expl.var.names # permet d'avoir les variables ordonnées dans le memes ordre
  
  # Calcul des courbes de réponse
  resp <- bm_PlotResponseCurves(bm.out = model_runs,
                                fixed.var = "mean", 
                                data_species = input_data@data.species # param permettant de donner les coordonnées des présences de l'espèce 
                                )$tab 
  # Courbes pour tous les md et run : difficilement lisible ! --> dans les lignes suivantes on refait un graph plus lisible
  # Faites la comparaison en enlevant le data_species !!
  
  # Renommage des colonnes  
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
  
  # A réfléchir : ajouter limites des données de présence ? --> bande noire sur l'axe des X du plot de biomod : permet de voir la densité des points de présences --> aide à l'interprétation --> permet d'expliquer que les courbes puissent être moins certaines quand il y a peu de points.
  
  
  png(paste0("graphiques/response_plot_", sp, ".png"), width = 550 * 4.2, height = 550 * 4.2, res = 300)
  
  print(p)
  
  dev.off()
}


# Interpretation des courbes pour Jaime.lannisterii (i=1)
# En bleu = moyenne de toutes les courbes des modèles = courbe approximative du md d'ensemble
# Biologiquement vraisemblable = pas de bimodalité qui suggèrerait une erreur sur nos modèles
# bio5 : sp n'a pas l'air très sensible : faible variation, et effectivement si on regarde son importance, elle est faible --> c'est peut-être une sp estivale qui est insensible aux températures hivernales, une sp qui tolère une grande gamme de temperature toute l'année

# Comparaison avec les "vraies courbes" dans /vs :
# bio5 : en général les vrai réponses à la température ressemble à ca : rapidement après l'optimum la temp devient mortelle et l'sp ne peut pas survivre. 

# bio5 : on ne retrouve pas la même courbe pcq on a pas échantillonné sur l'ensemble des température --> on observe d'ailleurs que derrière les md individuels ne sont pas du tout d'acc. 

# Interpretation des md indiv
# Dents de scie = pas ok = incohérente = on peut décider d'enlever ce md
# Bruit = normal
# Si courbes différentes selon les md, eg certaines logistique, certaines normales, certaines monotones positives --> se demander lesquels ont raison, pourquoi les modèles ne convergent pas ? 
# Très accidentée = modèles d'apprentissage automatique, si de grandes variations creux/vallées --> overfitting --> à éviter

# --> sélection manuelle des modèles 


# Avec Ned.starki (sp avec biais d'echant) on observe que les courbes individuelles sont très accidéntées et différentes.
