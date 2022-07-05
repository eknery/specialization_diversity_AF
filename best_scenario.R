setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")


best_geosse = read.table("geosse/geosse_best_fit_models.csv", sep=",", h = T)
best_geosse_time = read.table("geosse_time/geosse_time_best_fit_models.csv", sep=",", h = T)
best_quasse = read.table("quasse/quasse_best_fit_models.csv", sep=",", h = T)

best_scenario = c()
for (i in 1:100){
  aicc_values = c(best_geosse$aicc[i], best_geosse_time$aicc[i], best_quasse$aicc[i])
  names(aicc_values) = c("geosse", "geosse_time", "quasse")
  min_aicc = min(best_geosse$aicc[i], best_geosse_time$aicc[i], best_quasse$aicc[i])
  best_scenario = c( best_scenario, names(which(aicc_values == min_aicc)) )
}

best_scenario