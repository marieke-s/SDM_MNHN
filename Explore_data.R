# Load baseline data
a <- terra::rast("models/Jaime.lannisterii/proj_baseline/proj_baseline_Jaime.lannisterii.tif")

# Visu baseline projection across all models
plot(a)

# RUN1,2,3 correspond au run de validation croisées
# On observe que le md MARS prédit des valeurs + faible que les autre
# RF prédit des zones à forte probabilité plus petites que les autres md -> il a fait de l'overfitting (il faut regarder les courbes de réponses pour le voir) -> c'est un modèle dur à parametrer.
# GLM, GAM, GBoost sont bcp plus cohérents

# Visu modele d'ensemble : moyenne des proba de chaque md
plot(mean(a))

# Incertitude : standard deviation
sd_a <- app(a, sd) # Pour appliquer une fct sur toute les couches : app()
plot(sd_a)

# On observe une zone à forte variabilité entre les modèles dans l'Est : incertitude
# Pour comprendre ces zones de désaccord, il faut regarder les cartes individuellement etles courbes de réponses pour interpréter.
# On observe aussi des zones de faible incertitude : il s'agit de zones où les niches sont moins favorables (~ absence).

# Attention : ici on ne prend pas en compte l'incertitude lié à chaque modèle, ou l'incertitude au modèle climatique/ aux donnée --> il y'a de l'incertitude à différentes étapes mais c'est dur de répercuter l'incertitude à travers les étapes pour donner une incertitude totale.





