
iteration <- 100
method <- c("fofcd_fb_2_2", "glm")
  #"fofcd_fb", "fofcd_osb", "flars_basis", "gam", "glm"
simulation <- c("linear_1", "linear_2", "linear_5",
                "quadratic_5", "quadratic_10", "quadratic_20",
                "additive1_1", "additive1_2", "additive1_3",
                "additive2", "additive3",
                "interaction1", "interaction2")
for (imethod in 1:length(method)){
  filename <- c(paste0(rep("result_variable_selection_simu_", times = length(simulation)), 
                       simulation, "_", method[imethod], ".csv"))
  result <- data.frame(matrix(ncol = 3, nrow = length(filename)))
  for (idex in 1:length(filename)){
    data <- read.csv(filename[idex])
    coverage <- apply(data[, c(1,2)], 1, function(x) all(as.logical(x)))
    Smin <- apply(data[, c(1,2)], 1, function(x) ifelse(all(as.logical(x)==1), max(x), 0))
    result[idex, ] <- c(mean(coverage), mean(Smin[which(Smin!=0)]), mean(data[, ncol(data)]))
  }
  rownames(result) <- simulation
  colnames(result) <- c("Coverage", "Smin", "Time")
  write.csv(result, paste0("summary_variable_selection_simu_", method[imethod], ".csv"))
}          