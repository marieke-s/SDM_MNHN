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


# ----------- Species Range Change -----------
library(terra)

sp_list <- read.csv("data_cours/species_list.csv", sep = ";")


pdf("graphiques/src.pdf")
for (sp in sp_list$sp)
{
  src <- rast(paste("outputs/src_", sp, ".tif", sep = ""))
  for(j in names(src))
  {
    plot(src[[j]], legend = F, 
         main = paste0(sp, " range change\n", j), 
         col = c("#FF4100", "#3016B0", "#F2F2F2FF", "#2DD700"), xpd = NA, las = 1)
    legend("bottom", 
           fill = c("#F2F2F2FF", "#FF4100", "#3016B0", "#2DD700"), 
           legend = c("Unsuitable", "Lost", "Kept", "New"), ncol = 4, xpd = NA)
  }
}
dev.off()


pdf("graphiques/cont.pdf")
for(sp in sp_list$sp)
{
  probs <- rast(paste("outputs/em_cont_", sp, ".tif", sep = ""))
  plot(probs, main = paste(sp, names(probs)),
       col = rev(viridis::inferno(100)))
}
dev.off()

pdf("graphiques/binary.pdf")
for(sp in sp_list$sp)
{
  binary <- rast(paste("outputs/em_bin_", sp, ".tif", sep = ""))
  plot(binary, main = paste(sp, names(binary)))
}
dev.off()





# ---------- Figure standard deviations -----------
sp <- "Ned.starkii"
sd.sp <- rast(paste0("outputs/em_sd_", sp, ".tif"))
em.sp <- rast(paste0("outputs/em_cont_", sp, ".tif"))

cv.sp <- sd.sp/em.sp


cool = rainbow(50, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
lut  = c(rev(cool), rev(warm))

x11()
plot(cv.sp[["baseline"]], col = lut, main = "Coefficient of variation")
# Notez comme le CV semble pointer du doigt toutes les zones à faible proba
x11()
plot(sd.sp[["baseline"]], col = lut, main = "Standard deviation")



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


# ----------- Species Range Change -----------
library(terra)

sp_list <- read.csv("data_cours/species_list.csv", sep = ";")


pdf("graphiques/src.pdf")
for (sp in sp_list$sp)
{
  src <- rast(paste("outputs/src_", sp, ".tif", sep = ""))
  for(j in names(src))
  {
    plot(src[[j]], legend = F, 
         main = paste0(sp, " range change\n", j), 
         col = c("#FF4100", "#3016B0", "#F2F2F2FF", "#2DD700"), xpd = NA, las = 1)
    legend("bottom", 
           fill = c("#F2F2F2FF", "#FF4100", "#3016B0", "#2DD700"), 
           legend = c("Unsuitable", "Lost", "Kept", "New"), ncol = 4, xpd = NA)
  }
}
dev.off()


pdf("graphiques/cont.pdf")
for(sp in sp_list$sp)
{
  probs <- rast(paste("outputs/em_cont_", sp, ".tif", sep = ""))
  plot(probs, main = paste(sp, names(probs)),
       col = rev(viridis::inferno(100)))
}
dev.off()

pdf("graphiques/binary.pdf")
for(sp in sp_list$sp)
{
  binary <- rast(paste("outputs/em_bin_", sp, ".tif", sep = ""))
  plot(binary, main = paste(sp, names(binary)))
}
dev.off()





# ---------- Figure standard deviations -----------
sp <- "Ned.starkii"
sd.sp <- rast(paste0("outputs/em_sd_", sp, ".tif"))
em.sp <- rast(paste0("outputs/em_cont_", sp, ".tif"))

cv.sp <- sd.sp/em.sp


cool = rainbow(50, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
lut  = c(rev(cool), rev(warm))

x11()
plot(cv.sp[["baseline"]], col = lut, main = "Coefficient of variation")
# Notez comme le CV semble pointer du doigt toutes les zones à faible proba
x11()
plot(sd.sp[["baseline"]], col = lut, main = "Standard deviation")



>>>>>>> 5768564 (2nd commit)
