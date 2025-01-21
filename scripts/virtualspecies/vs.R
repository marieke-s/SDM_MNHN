<<<<<<< HEAD
library(virtualspecies) 
library(geometry)
library(rgl)

# Function to grab midpoints of intervals
midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
  upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
  return(round(lower+(upper-lower)/2, dp))
}





baseline <- stack("data/baseline.tif")
a <- formatFunctions(baseline[[c("bio5", "bio6", "forests")]],
                     bio5 = c(fun = "betaFun", p1 = 9, p2 = 34, alpha = 1, gamma = 0.15),
                     bio6 = c(fun = "custnorm", mean = -5, diff = 10, prob = .99),
                     forests = c(fun = "logisticFun", beta = 40, alpha = -6))

sp.virtuelle <- generateSpFromFun(baseline[[c("bio5", "bio6", "forests")]],
                                    parameters = a, plot = TRUE)

sp.virtuelle <- convertToPA(sp.virtuelle, prob.method = "linear")

png("./data_cours/vs/vs_relations.png", w = 900, h = 900)
plotResponse(sp.virtuelle)
dev.off()

png("./data_cours/vs/vs_distribution.png", w = 1200, h = 700)
plot(sp.virtuelle)
dev.off()

saveRDS(sp.virtuelle,
        "./data_cours/virtualspecies.rds")


# 1. 2000 P/A, random sampling
windows(6000, 4000)
Jaime.lannisterii <- sampleOccurrences(sp.virtuelle, 
                                      n = 2000, type = "presence-absence",
                                      extract.probability = TRUE,
                                      replacement = TRUE)$sample.points


# 2. 100 PA, random sampling
Tyrion.lannisterii <- sampleOccurrences(sp.virtuelle, 
                                        n = 100, type = "presence-absence",
                                        extract.probability = TRUE,
                                        replacement = TRUE)$sample.points

# 3. 1000 P, random sampling
Jon.snowii <- sampleOccurrences(sp.virtuelle, 
                                      n = 1000,
                                      extract.probability = TRUE,
                                      replacement = TRUE)$sample.points

sampling.bias <- readRDS("data_cours/vs/distance_roads.RDS")
sampling.bias <- resample(sampling.bias, baseline)
sampling.bias <- sampling.bias@data@max - sampling.bias

# 4. 1000 P, sampling bias
Ned.starkii <- sampleOccurrences(sp.virtuelle, 
                                 n = 100,
                                 bias = "manual",
                                 weights = sampling.bias,
                                 bias.area = "France",
                                 detection.probability = 0.5,
                                 extract.probability = TRUE,
                                 replacement = TRUE)$sample.points
Ned.starkii.obs <- Ned.starkii[-which(is.na(Ned.starkii$Observed)), ]
any(duplicated(Ned.starkii.obs))

Hermione.grangerii <- sampleOccurrences(sp.virtuelle, 
                                    n = 50,
                                    extract.probability = TRUE,
                                    replacement = FALSE)$sample.points


write.table(Hermione.grangerii[, c("x", "y", "Observed")],
            "./data/Hermione.grangerii.csv", sep  = ";")

write.table(Jaime.lannisterii[, c("x", "y", "Observed", "true.probability")],
            "./data_cours/Jaime.lannisterii.csv", sep  = ";")
write.table(Tyrion.lannisterii[, c("x", "y", "Observed", "true.probability")],
            "./data_cours/Tyrion.lannisterii.csv", sep  = ";")
write.table(Jon.snowii[, c("x", "y", "Observed", "true.probability")],
            "./data_cours/Jon.snowii.csv", sep  = ";")
write.table(Ned.starkii.obs[, c("x", "y", "Observed", "true.probability")],
            "./data_cours/Ned.starkii.csv", sep  = ";")

sp.list <- data.frame(sp = c("Jaime.lannisterii", 
                             "Tyrion.lannisterii",
                             "Jon.snowii",
                             "Ned.starkii"),
                      occ = c(rep("PA", 2), rep("PO", 2)),
                      env.data.type = c(rep("Spatial", 4)),
                      pa.generation = c(rep("no pseudoabs", 2),
                                        rep("biomod", 2)))
write.table(sp.list, 
            "./data_cours/species_list.csv", sep = ";")


future_data <- stack("./data/ipsl_ssp245_2070.tif")
sp.virtuelle.2.6 <- generateSpFromFun(future_data[[c("bio5", "bio6", "forests")]],
                                  parameters = a, plot = TRUE)
sp.virtuelle.2.6 <- convertToPA(sp.virtuelle.2.6, prob.method = "linear")

png("./data_cours/vs/vs_distribution_futur2.6.png", w = 1200, h = 700)
plot(sp.virtuelle.2.6)
dev.off()

src.2.6 <- sp.virtuelle.2.6$pa.raster - 2 * sp.virtuelle$pa.raster
png("./data_cours/vs/vs_src_2.6.png", w = 900, h = 900)
plot(src.2.6, legend = F, 
     main = "Species range change 2.6 2050", 
     col = c("#FF4100", "#3016B0", "#F2F2F2FF", "#2DD700"), xpd = NA, las = 1)
legend("bottom", 
       fill = c("#F2F2F2FF", "#FF4100", "#3016B0", "#2DD700"), 
       legend = c("Unsuitable", "Lost", "Kept", "New"), ncol = 4, xpd = NA)
dev.off()


future_data <- stack("./data/ipsl_ssp585_2070.tif")
sp.virtuelle.8.5 <- generateSpFromFun(future_data[[c("bio5", "bio6", "forests")]],
                                  parameters = a, plot = TRUE)
sp.virtuelle.8.5 <- convertToPA(sp.virtuelle.8.5, prob.method = "linear")

png("./data_cours/vs/vs_distribution_futur8.5.png", w = 1200, h = 700)
plot(sp.virtuelle.8.5)
dev.off()

src.8.5 <- sp.virtuelle.8.5$pa.raster - 2 * sp.virtuelle$pa.raster
png("./data_cours/vs/vs_src_8.5.png", w = 900, h = 900)
plot(src.8.5, legend = F, 
     main = "Species range change 8.5 2050", 
     col = c("#FF4100", "#3016B0", "#F2F2F2FF", "#2DD700"), xpd = NA, las = 1)
legend("bottom", 
       fill = c("#F2F2F2FF", "#FF4100", "#3016B0", "#2DD700"), 
       legend = c("Unsuitable", "Lost", "Kept", "New"), ncol = 4, xpd = NA)
dev.off()


# # Environmental conditions
# combinations <- getValues(baseline[[c("bio5", "bio6", "forests")]])
# all.duplicated.conditions <- which(is.na(combinations[, 1]))
# all.xy <- xyFromCell(baseline, 1:ncell(baseline))
# all.xy <- all.xy[-all.duplicated.conditions, ]
# combinations <- combinations[-all.duplicated.conditions, ]
# 
# min.max <- data.frame(min = apply(combinations, 2, min),
#                       max = apply(combinations, 2, max)) 
# 
# 
# # Creation of the intervals for each variable
# # We use intervals specific to each variable
# # Should consider adding max and min independently for some variables when they have
# # a specific meaning
# # Maybe also think about log transforming variables when they are skewed (e.g. precipitation)
# seqs <- list(bio5 = seq(min.max["bio5", "min"], min.max["bio5", "max"], length = 50),
#              bio6 = seq(min.max["bio6", "min"], min.max["bio6", "max"], length = 50),
#              forests = seq(min.max["forests", "min"], min.max["forests", "max"], length = 11))
# 
# 
# # All possible combination in the environment
# comb.cat <- sapply(colnames(combinations), function(x, combs, seqs.)
#   cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
#   combs = combinations, seqs. = seqs)      
# message("Total number of cells with environmental conditions in the geographical space: ", nrow(combinations),
#         "\nNumber of duplicated conditions: ", length(which(duplicated(comb.cat))),
#         "\nNumber of unique cells (environmental space): ", nrow(comb.cat[-which(duplicated(comb.cat)), ]))
# possible.combs <- lapply(seqs, function(x)
#   data.frame(interval = cut(x, x, right = FALSE, include.lowest = TRUE),
#              mids = midpoints(cut(x, x, right = FALSE, include.lowest = TRUE))))
# 
# all.env.pixels <- sapply(colnames(comb.cat), function(x, int.to.replace, replacing.values)
# {
#   replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
# }, int.to.replace = comb.cat, replacing.values = possible.combs)
# 
# duplicated.cells <- which(duplicated(all.env.pixels))
# if(length(duplicated.cells))
# {
#   all.env.pixels <- all.env.pixels[-duplicated.cells, ]
#   all.xy <- all.xy[-duplicated.cells, ]
# }
# plot3d(all.env.pixels) # First 3 axes
# 
# plot(sp.virtuelle$suitab.raster)
# plot(baseline[["bio5"]])
# points(all.xy, cex = .1)
# 
# 
# 
# # 9. Kit filtrage environnemental + pseudoabsences
# env.conditionskit <- extract(baseline[[c("bio5", "bio6","forests")]],
#                              Kit.harringtonii[, c("x", "y")])
# 
# env.combkit <- sapply(colnames(env.conditionskit), function(x, combs, seqs.)
#   cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
#   combs = env.conditionskit, seqs. = seqs)
# 
# 
# message("Total number of cells with environmental conditions in the geographical space: ", nrow(env.combkit),
#         "\nNumber of duplicated conditions: ", length(which(duplicated(env.combkit))),
#         "\nNumber of unique cells (environmental space): ", nrow(env.combkit[-which(duplicated(env.combkit)), ]))
# 
# duplicated.cells <- which(duplicated(env.combkit))
# if(length(duplicated.cells))
# {
#   env.combkit <- cbind(env.combkit[-duplicated.cells, ])
#   Bruce.willis <- cbind(Kit.harringtonii[-duplicated.cells, c("x", "y", "Observed", "true.probability")],
#                         sapply(colnames(env.combkit), function(x, int.to.replace, replacing.values)
#                         {
#                           replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
#                         }, int.to.replace = env.combkit, replacing.values = possible.combs))
# } else
# {
#   Bruce.willis <- cbind(Kit.harringtonii[, c("x", "y", "Observed", "true.probability")],
#                     sapply(colnames(env.combkit), function(x, int.to.replace, replacing.values)
#                     {
#                       replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
#                     }, int.to.replace = env.combkit, replacing.values = possible.combs))
# }
# 
# 
# # Calculating the convex hull of environmental conditions
# Bruce.willis.matrix <- as.matrix(Bruce.willis[, c("bio5", "bio6", "forests")])
# Bruce.willis.convhull <- convhulln(Bruce.willis.matrix) 
# 
# # Checking environmental conditions that are inside the convex hull vs. outside
# Bruce.willis.inhull <- inhulln(Bruce.willis.convhull, 
#                      all.env.pixels)
# 
# 
# 
# plot3d(all.env.pixels[!Bruce.willis.inhull, c("bio5", "bio6", "forests")],
#        col = "green", size = 5) # First 3 axes
# plot3d(all.env.pixels[Bruce.willis.inhull, c("bio5", "bio6", "forests")],
#        col = "red", add = T, size = 10) # First 3 axes
# plot3d(Bruce.willis.matrix[, c("bio5", "bio6", "forests")], add = TRUE, size = 15)
# # NOTE: some points can be both IN and out of the convex hull
# # because of the fourth dimension
# # They are green+red in the graph
# 
# 
# # Sampling pseudo-absences outside the convex hull
# Bruce.willis.pseudoabs <- sample(which(!Bruce.willis.inhull), 
#                                  size = nrow(Bruce.willis),
#                                  replace = FALSE)
# Bruce.willis <- rbind.data.frame(Bruce.willis,
#                                  data.frame(all.xy[Bruce.willis.pseudoabs, ],
#                                             Observed = 0,
#                                             true.probability = extract(sp.virtuelle$probability.of.occurrence,
#                                                                        all.xy[Bruce.willis.pseudoabs, ]),
#                                             all.env.pixels[Bruce.willis.pseudoabs, ]))
# 
# 
# plot(sp.virtuelle$pa.raster)
# points(Bruce.willis[, c("x", "y")],
#        pch = c(1, 16)[as.factor(Bruce.willis$Observed)],
#        cex=.5)
# 
# # 10. Sean filtrage environnemental + pseudoabsences
# Sean.beanus.obs <- Sean.beanus[-which(is.na(Sean.beanus$Observed)), ]
# env.conditionsSean <- extract(baseline[[c("bio5", "bio6", "forests")]],
#                               Sean.beanus.obs[, c("x", "y")])
# 
# env.combSean <- sapply(colnames(env.conditionsSean), function(x, combs, seqs.)
#   cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
#   combs = env.conditionsSean, seqs. = seqs)
# 
# 
# message("Total number of cells with environmental conditions in the geographical space: ", nrow(env.combSean),
#         "\nNumber of duplicated conditions: ", length(which(duplicated(env.combSean))),
#         "\nNumber of unique cells (environmental space): ", nrow(env.combSean[-which(duplicated(env.combSean)), ]))
# 
# duplicated.cells <- which(duplicated(env.combSean))
# if(length(duplicated.cells))
# {
#   env.combSean <- cbind(env.combSean[-duplicated.cells, ])
#   Iain.glainus <- cbind(Sean.beanus.obs[-duplicated.cells, c("x", "y", "Observed", "true.probability")],
#                 sapply(colnames(env.combSean), function(x, int.to.replace, replacing.values)
#   {
#     replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
#   }, int.to.replace = env.combSean, replacing.values = possible.combs))
# } else
# {
#   Iain.glainus <- cbind(Sean.beanus.obs[, c("x", "y", "Observed", "true.probability")],
#                     sapply(colnames(env.combSean), function(x, int.to.replace, replacing.values)
#   {
#     replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
#   }, int.to.replace = env.combSean, replacing.values = possible.combs))
# }
# 
# 
# # Calculating the convex hull of environmental conditions
# Iain.glainus.matrix <- as.matrix(Iain.glainus[, c("bio5", "bio6", "forests")])
# Iain.glainus.convhull <- convhulln(Iain.glainus.matrix) 
# 
# # Checking environmental conditions that are inside the convex hull vs. outside
# Iain.glainus.inhull <- inhulln(Iain.glainus.convhull, 
#                      all.env.pixels)
# 
# 
# plot3d(all.env.pixels[!Iain.glainus.inhull, c("bio5", "bio6", "forests")],
#        col = "green")
# plot3d(all.env.pixels[Iain.glainus.inhull, c("bio5", "bio6", "forests")],
#        col = "red", size = 10, add = T)
# plot3d(Iain.glainus[, c("bio5", "bio6", "forests")], size = 15, add = T)
# # NOTE: some points can be both IN and out of the convex hull
# # because of the fourth dimension
# # They are green+red in the graph
# 
# 
# # Sampling pseudo-absences outside the convex hull
# Iain.glainus.pseudoabs <- sample(which(!Iain.glainus.inhull), 
#                        size = nrow(Iain.glainus),
#                        replace = FALSE)
# Iain.glainus <- rbind.data.frame(Iain.glainus,
#                              data.frame(all.xy[Iain.glainus.pseudoabs, ],
#                                         Observed = 0,
#                                         true.probability = extract(sp.virtuelle$probability.of.occurrence,
#                                                                    all.xy[Iain.glainus.pseudoabs, ]),
#                                         all.env.pixels[Iain.glainus.pseudoabs, ]))
# 
# plot(sp.virtuelle$pa.raster)
# points(Iain.glainus[, c("x", "y")],
#        pch = c(1, 16)[as.factor(Iain.glainus$Observed)],
#        cex=.5)
# 
# write.table(Bruce.willis,
#             "./data_cours/Bruce.willis.csv", sep  = ";")
# write.table(Iain.glainus,
#             "./data_cours/Iain.glainus.csv", sep  = ";")
# 
# saveRDS(Bruce.willis, "./data_cours/occurrences_Bruce.willis")
# saveRDS(Iain.glainus, "./data_cours/occurrences_Iain.glainus")

=======
library(virtualspecies) 
library(geometry)
library(rgl)

# Function to grab midpoints of intervals
midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
  upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
  return(round(lower+(upper-lower)/2, dp))
}





baseline <- stack("data/baseline.tif")
a <- formatFunctions(baseline[[c("bio5", "bio6", "forests")]],
                     bio5 = c(fun = "betaFun", p1 = 9, p2 = 34, alpha = 1, gamma = 0.15),
                     bio6 = c(fun = "custnorm", mean = -5, diff = 10, prob = .99),
                     forests = c(fun = "logisticFun", beta = 40, alpha = -6))

sp.virtuelle <- generateSpFromFun(baseline[[c("bio5", "bio6", "forests")]],
                                    parameters = a, plot = TRUE)

sp.virtuelle <- convertToPA(sp.virtuelle, prob.method = "linear")

png("./data_cours/vs/vs_relations.png", w = 900, h = 900)
plotResponse(sp.virtuelle)
dev.off()

png("./data_cours/vs/vs_distribution.png", w = 1200, h = 700)
plot(sp.virtuelle)
dev.off()

saveRDS(sp.virtuelle,
        "./data_cours/virtualspecies.rds")


# 1. 2000 P/A, random sampling
windows(6000, 4000)
Jaime.lannisterii <- sampleOccurrences(sp.virtuelle, 
                                      n = 2000, type = "presence-absence",
                                      extract.probability = TRUE,
                                      replacement = TRUE)$sample.points


# 2. 100 PA, random sampling
Tyrion.lannisterii <- sampleOccurrences(sp.virtuelle, 
                                        n = 100, type = "presence-absence",
                                        extract.probability = TRUE,
                                        replacement = TRUE)$sample.points

# 3. 1000 P, random sampling
Jon.snowii <- sampleOccurrences(sp.virtuelle, 
                                      n = 1000,
                                      extract.probability = TRUE,
                                      replacement = TRUE)$sample.points

sampling.bias <- readRDS("data_cours/vs/distance_roads.RDS")
sampling.bias <- resample(sampling.bias, baseline)
sampling.bias <- sampling.bias@data@max - sampling.bias

# 4. 1000 P, sampling bias
Ned.starkii <- sampleOccurrences(sp.virtuelle, 
                                 n = 100,
                                 bias = "manual",
                                 weights = sampling.bias,
                                 bias.area = "France",
                                 detection.probability = 0.5,
                                 extract.probability = TRUE,
                                 replacement = TRUE)$sample.points
Ned.starkii.obs <- Ned.starkii[-which(is.na(Ned.starkii$Observed)), ]
any(duplicated(Ned.starkii.obs))

Hermione.grangerii <- sampleOccurrences(sp.virtuelle, 
                                    n = 50,
                                    extract.probability = TRUE,
                                    replacement = FALSE)$sample.points


write.table(Hermione.grangerii[, c("x", "y", "Observed")],
            "./data/Hermione.grangerii.csv", sep  = ";")

write.table(Jaime.lannisterii[, c("x", "y", "Observed", "true.probability")],
            "./data_cours/Jaime.lannisterii.csv", sep  = ";")
write.table(Tyrion.lannisterii[, c("x", "y", "Observed", "true.probability")],
            "./data_cours/Tyrion.lannisterii.csv", sep  = ";")
write.table(Jon.snowii[, c("x", "y", "Observed", "true.probability")],
            "./data_cours/Jon.snowii.csv", sep  = ";")
write.table(Ned.starkii.obs[, c("x", "y", "Observed", "true.probability")],
            "./data_cours/Ned.starkii.csv", sep  = ";")

sp.list <- data.frame(sp = c("Jaime.lannisterii", 
                             "Tyrion.lannisterii",
                             "Jon.snowii",
                             "Ned.starkii"),
                      occ = c(rep("PA", 2), rep("PO", 2)),
                      env.data.type = c(rep("Spatial", 4)),
                      pa.generation = c(rep("no pseudoabs", 2),
                                        rep("biomod", 2)))
write.table(sp.list, 
            "./data_cours/species_list.csv", sep = ";")


future_data <- stack("./data/ipsl_ssp245_2070.tif")
sp.virtuelle.2.6 <- generateSpFromFun(future_data[[c("bio5", "bio6", "forests")]],
                                  parameters = a, plot = TRUE)
sp.virtuelle.2.6 <- convertToPA(sp.virtuelle.2.6, prob.method = "linear")

png("./data_cours/vs/vs_distribution_futur2.6.png", w = 1200, h = 700)
plot(sp.virtuelle.2.6)
dev.off()

src.2.6 <- sp.virtuelle.2.6$pa.raster - 2 * sp.virtuelle$pa.raster
png("./data_cours/vs/vs_src_2.6.png", w = 900, h = 900)
plot(src.2.6, legend = F, 
     main = "Species range change 2.6 2050", 
     col = c("#FF4100", "#3016B0", "#F2F2F2FF", "#2DD700"), xpd = NA, las = 1)
legend("bottom", 
       fill = c("#F2F2F2FF", "#FF4100", "#3016B0", "#2DD700"), 
       legend = c("Unsuitable", "Lost", "Kept", "New"), ncol = 4, xpd = NA)
dev.off()


future_data <- stack("./data/ipsl_ssp585_2070.tif")
sp.virtuelle.8.5 <- generateSpFromFun(future_data[[c("bio5", "bio6", "forests")]],
                                  parameters = a, plot = TRUE)
sp.virtuelle.8.5 <- convertToPA(sp.virtuelle.8.5, prob.method = "linear")

png("./data_cours/vs/vs_distribution_futur8.5.png", w = 1200, h = 700)
plot(sp.virtuelle.8.5)
dev.off()

src.8.5 <- sp.virtuelle.8.5$pa.raster - 2 * sp.virtuelle$pa.raster
png("./data_cours/vs/vs_src_8.5.png", w = 900, h = 900)
plot(src.8.5, legend = F, 
     main = "Species range change 8.5 2050", 
     col = c("#FF4100", "#3016B0", "#F2F2F2FF", "#2DD700"), xpd = NA, las = 1)
legend("bottom", 
       fill = c("#F2F2F2FF", "#FF4100", "#3016B0", "#2DD700"), 
       legend = c("Unsuitable", "Lost", "Kept", "New"), ncol = 4, xpd = NA)
dev.off()


# # Environmental conditions
# combinations <- getValues(baseline[[c("bio5", "bio6", "forests")]])
# all.duplicated.conditions <- which(is.na(combinations[, 1]))
# all.xy <- xyFromCell(baseline, 1:ncell(baseline))
# all.xy <- all.xy[-all.duplicated.conditions, ]
# combinations <- combinations[-all.duplicated.conditions, ]
# 
# min.max <- data.frame(min = apply(combinations, 2, min),
#                       max = apply(combinations, 2, max)) 
# 
# 
# # Creation of the intervals for each variable
# # We use intervals specific to each variable
# # Should consider adding max and min independently for some variables when they have
# # a specific meaning
# # Maybe also think about log transforming variables when they are skewed (e.g. precipitation)
# seqs <- list(bio5 = seq(min.max["bio5", "min"], min.max["bio5", "max"], length = 50),
#              bio6 = seq(min.max["bio6", "min"], min.max["bio6", "max"], length = 50),
#              forests = seq(min.max["forests", "min"], min.max["forests", "max"], length = 11))
# 
# 
# # All possible combination in the environment
# comb.cat <- sapply(colnames(combinations), function(x, combs, seqs.)
#   cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
#   combs = combinations, seqs. = seqs)      
# message("Total number of cells with environmental conditions in the geographical space: ", nrow(combinations),
#         "\nNumber of duplicated conditions: ", length(which(duplicated(comb.cat))),
#         "\nNumber of unique cells (environmental space): ", nrow(comb.cat[-which(duplicated(comb.cat)), ]))
# possible.combs <- lapply(seqs, function(x)
#   data.frame(interval = cut(x, x, right = FALSE, include.lowest = TRUE),
#              mids = midpoints(cut(x, x, right = FALSE, include.lowest = TRUE))))
# 
# all.env.pixels <- sapply(colnames(comb.cat), function(x, int.to.replace, replacing.values)
# {
#   replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
# }, int.to.replace = comb.cat, replacing.values = possible.combs)
# 
# duplicated.cells <- which(duplicated(all.env.pixels))
# if(length(duplicated.cells))
# {
#   all.env.pixels <- all.env.pixels[-duplicated.cells, ]
#   all.xy <- all.xy[-duplicated.cells, ]
# }
# plot3d(all.env.pixels) # First 3 axes
# 
# plot(sp.virtuelle$suitab.raster)
# plot(baseline[["bio5"]])
# points(all.xy, cex = .1)
# 
# 
# 
# # 9. Kit filtrage environnemental + pseudoabsences
# env.conditionskit <- extract(baseline[[c("bio5", "bio6","forests")]],
#                              Kit.harringtonii[, c("x", "y")])
# 
# env.combkit <- sapply(colnames(env.conditionskit), function(x, combs, seqs.)
#   cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
#   combs = env.conditionskit, seqs. = seqs)
# 
# 
# message("Total number of cells with environmental conditions in the geographical space: ", nrow(env.combkit),
#         "\nNumber of duplicated conditions: ", length(which(duplicated(env.combkit))),
#         "\nNumber of unique cells (environmental space): ", nrow(env.combkit[-which(duplicated(env.combkit)), ]))
# 
# duplicated.cells <- which(duplicated(env.combkit))
# if(length(duplicated.cells))
# {
#   env.combkit <- cbind(env.combkit[-duplicated.cells, ])
#   Bruce.willis <- cbind(Kit.harringtonii[-duplicated.cells, c("x", "y", "Observed", "true.probability")],
#                         sapply(colnames(env.combkit), function(x, int.to.replace, replacing.values)
#                         {
#                           replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
#                         }, int.to.replace = env.combkit, replacing.values = possible.combs))
# } else
# {
#   Bruce.willis <- cbind(Kit.harringtonii[, c("x", "y", "Observed", "true.probability")],
#                     sapply(colnames(env.combkit), function(x, int.to.replace, replacing.values)
#                     {
#                       replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
#                     }, int.to.replace = env.combkit, replacing.values = possible.combs))
# }
# 
# 
# # Calculating the convex hull of environmental conditions
# Bruce.willis.matrix <- as.matrix(Bruce.willis[, c("bio5", "bio6", "forests")])
# Bruce.willis.convhull <- convhulln(Bruce.willis.matrix) 
# 
# # Checking environmental conditions that are inside the convex hull vs. outside
# Bruce.willis.inhull <- inhulln(Bruce.willis.convhull, 
#                      all.env.pixels)
# 
# 
# 
# plot3d(all.env.pixels[!Bruce.willis.inhull, c("bio5", "bio6", "forests")],
#        col = "green", size = 5) # First 3 axes
# plot3d(all.env.pixels[Bruce.willis.inhull, c("bio5", "bio6", "forests")],
#        col = "red", add = T, size = 10) # First 3 axes
# plot3d(Bruce.willis.matrix[, c("bio5", "bio6", "forests")], add = TRUE, size = 15)
# # NOTE: some points can be both IN and out of the convex hull
# # because of the fourth dimension
# # They are green+red in the graph
# 
# 
# # Sampling pseudo-absences outside the convex hull
# Bruce.willis.pseudoabs <- sample(which(!Bruce.willis.inhull), 
#                                  size = nrow(Bruce.willis),
#                                  replace = FALSE)
# Bruce.willis <- rbind.data.frame(Bruce.willis,
#                                  data.frame(all.xy[Bruce.willis.pseudoabs, ],
#                                             Observed = 0,
#                                             true.probability = extract(sp.virtuelle$probability.of.occurrence,
#                                                                        all.xy[Bruce.willis.pseudoabs, ]),
#                                             all.env.pixels[Bruce.willis.pseudoabs, ]))
# 
# 
# plot(sp.virtuelle$pa.raster)
# points(Bruce.willis[, c("x", "y")],
#        pch = c(1, 16)[as.factor(Bruce.willis$Observed)],
#        cex=.5)
# 
# # 10. Sean filtrage environnemental + pseudoabsences
# Sean.beanus.obs <- Sean.beanus[-which(is.na(Sean.beanus$Observed)), ]
# env.conditionsSean <- extract(baseline[[c("bio5", "bio6", "forests")]],
#                               Sean.beanus.obs[, c("x", "y")])
# 
# env.combSean <- sapply(colnames(env.conditionsSean), function(x, combs, seqs.)
#   cut(combs[, x], breaks = seqs.[[x]], include.lowest = TRUE, right = FALSE),
#   combs = env.conditionsSean, seqs. = seqs)
# 
# 
# message("Total number of cells with environmental conditions in the geographical space: ", nrow(env.combSean),
#         "\nNumber of duplicated conditions: ", length(which(duplicated(env.combSean))),
#         "\nNumber of unique cells (environmental space): ", nrow(env.combSean[-which(duplicated(env.combSean)), ]))
# 
# duplicated.cells <- which(duplicated(env.combSean))
# if(length(duplicated.cells))
# {
#   env.combSean <- cbind(env.combSean[-duplicated.cells, ])
#   Iain.glainus <- cbind(Sean.beanus.obs[-duplicated.cells, c("x", "y", "Observed", "true.probability")],
#                 sapply(colnames(env.combSean), function(x, int.to.replace, replacing.values)
#   {
#     replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
#   }, int.to.replace = env.combSean, replacing.values = possible.combs))
# } else
# {
#   Iain.glainus <- cbind(Sean.beanus.obs[, c("x", "y", "Observed", "true.probability")],
#                     sapply(colnames(env.combSean), function(x, int.to.replace, replacing.values)
#   {
#     replacing.values[[x]]$mids[match(int.to.replace[, x], replacing.values[[x]]$interval)]
#   }, int.to.replace = env.combSean, replacing.values = possible.combs))
# }
# 
# 
# # Calculating the convex hull of environmental conditions
# Iain.glainus.matrix <- as.matrix(Iain.glainus[, c("bio5", "bio6", "forests")])
# Iain.glainus.convhull <- convhulln(Iain.glainus.matrix) 
# 
# # Checking environmental conditions that are inside the convex hull vs. outside
# Iain.glainus.inhull <- inhulln(Iain.glainus.convhull, 
#                      all.env.pixels)
# 
# 
# plot3d(all.env.pixels[!Iain.glainus.inhull, c("bio5", "bio6", "forests")],
#        col = "green")
# plot3d(all.env.pixels[Iain.glainus.inhull, c("bio5", "bio6", "forests")],
#        col = "red", size = 10, add = T)
# plot3d(Iain.glainus[, c("bio5", "bio6", "forests")], size = 15, add = T)
# # NOTE: some points can be both IN and out of the convex hull
# # because of the fourth dimension
# # They are green+red in the graph
# 
# 
# # Sampling pseudo-absences outside the convex hull
# Iain.glainus.pseudoabs <- sample(which(!Iain.glainus.inhull), 
#                        size = nrow(Iain.glainus),
#                        replace = FALSE)
# Iain.glainus <- rbind.data.frame(Iain.glainus,
#                              data.frame(all.xy[Iain.glainus.pseudoabs, ],
#                                         Observed = 0,
#                                         true.probability = extract(sp.virtuelle$probability.of.occurrence,
#                                                                    all.xy[Iain.glainus.pseudoabs, ]),
#                                         all.env.pixels[Iain.glainus.pseudoabs, ]))
# 
# plot(sp.virtuelle$pa.raster)
# points(Iain.glainus[, c("x", "y")],
#        pch = c(1, 16)[as.factor(Iain.glainus$Observed)],
#        cex=.5)
# 
# write.table(Bruce.willis,
#             "./data_cours/Bruce.willis.csv", sep  = ";")
# write.table(Iain.glainus,
#             "./data_cours/Iain.glainus.csv", sep  = ";")
# 
# saveRDS(Bruce.willis, "./data_cours/occurrences_Bruce.willis")
# saveRDS(Iain.glainus, "./data_cours/occurrences_Iain.glainus")

>>>>>>> 5768564 (2nd commit)
