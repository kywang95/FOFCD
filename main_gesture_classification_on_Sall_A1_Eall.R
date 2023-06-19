library(FuncNN)
library(fda)
library(ggpubr)
library(keras)
library(RColorBrewer)
library(flars)
source("function_FOFCD.R")

nsubject <- 27
nsensor <- 10+22
# nclass_1 <- 12
# nclass_2 <- 17
# nclass_3 <- 23
# nclass <- 12 # 12 # 17 # 8 # 9 # 23 # 12+17+23
nrep <- 10
simulation <- c("E1", "E2", "E3", "Eall")
classes <- c(12, 17, 23, 12+17+23)

for (i_simu in 1:length(simulation)){
  cat("simulation_", simulation[i_simu], "\n")
  # i_simu <- 1
  nclass <- classes[i_simu]
  train_n <- nsubject*nclass*nrep*0.7
  test_n <- nsubject*nclass*nrep*0.3
  
  ydata <- read.csv(paste0("Sall_A1_", simulation[i_simu], "_class_of_posture.csv"))
  y <- ydata$class
  train_y <- y[rep(
    rep((1:nclass), times = nrep*0.7) + 
      rep((c(1,3,4,6,8,9,10)-1)*nclass, each = nclass), 
    times = nsubject) + 
      rep(((1:nsubject) - 1)*nclass*nrep, each = train_n/nsubject)]
  test_y <- y[rep(
    rep((1:nclass), times = nrep*0.3) + 
      rep((c(2,5,7)-1)*nclass, each = nclass), 
    times = nsubject) + 
      rep(((1:nsubject) - 1)*nclass*nrep, each = test_n/nsubject)]
  
  train_signal <- list()
  test_signal <- list()
  for(isensor in 1:nsensor){
    signal <- as.matrix(read.csv(paste0("Sall_A1_", simulation[i_simu], "_signal_of_sensor_", isensor, ".csv")))[,1:500]
    row.na <- which(is.na(signal[, ncol(signal)]))
    if (length(row.na)){
      for (irow in 1:length(row.na)){
        col.na <- which(is.na(signal[row.na[irow], ]))
        signal[row.na[irow], col.na] <- 0
      }
    }
    train_signal[[isensor]] <- signal[rep(
      rep((1:nclass), times = nrep*0.7) + 
        rep((c(1,3,4,6,8,9,10)-1)*nclass, each = nclass), 
      times = nsubject) + 
        rep(((1:nsubject) - 1)*nclass*nrep, each = train_n/nsubject), ]
    test_signal[[isensor]] <- signal[rep(
      rep((1:nclass), times = nrep*0.3) + 
        rep((c(2,5,7)-1)*nclass, each = nclass), 
      times = nsubject) + 
        rep(((1:nsubject) - 1)*nclass*nrep, each = test_n/nsubject), ]
  }
  
  ntime <- ncol(train_signal[[1]])
  TRange=c(0, 1)
  tmax = max(TRange)
  tmin = min(TRange)
  Time <- seq(from = tmin, to = tmax, length.out = ntime)
  
  # define the fourier basis 
  nbasis = 65
  temp_basis = create.fourier.basis(TRange, nbasis)
  timepts <- Time
  # convert the functional predictor in training data into fda objects
  train_fda = array(dim = c(nbasis, train_n, nsensor))
  test_fda = array(dim = c(nbasis, test_n, nsensor))
  train_data <- list()
  test_data <- list()
  for(isensor in 1:nsensor){
    train_signal_data =  Data2fd(timepts, t(train_signal[[isensor]]), temp_basis)
    train_data[[isensor]] <- t(eval.fd(Time, train_signal_data))
    func_cov = train_signal_data$coefs
    train_fda[,,isensor] <- func_cov
    test_signal_data =  Data2fd(timepts, t(test_signal[[isensor]]), temp_basis)
    test_data[[isensor]] <- t(eval.fd(Time, test_signal_data))
    func_cov = test_signal_data$coefs
    test_fda[,,isensor] <- func_cov
  }
  
  data <- list(y = matrix(train_y, ncol = 1), x = train_data, Time = Time)
  
  ## FOFCD
  method = 'fourier_basis'; control = list(nbasis = 4, norder = 4, t = Time)
  fofcd.fb <- fofcd(data$y, data$x, numCores = 1, method = method, control = control)
  # print(fofcd.fb)
  feature.fb <- fofcd.fb$feature$col_remain
  print(feature.fb)
  print(fofcd.fb$time)
  write.csv(feature.fb, paste0("feature_fofcd_", simulation[i_simu], ".csv"))
  
  # ## fLARS
  # flars.basis <- flars::flars(data$x, data$y, method='basis', 
  #                             control = list(t = timepts), 
  #                             max_selection = nsensor,
  #                             normalize='norm'
  #                             # , lasso=TRUE
  # )
  # feature.fLARS <- flars.basis$index
  # print(feature.fLARS)
  # if (i_simu == 1){
  #   feature.fLARS <- c(16, 15, 13, 25, 14, 21,32, 12, 7)
  # }else
  #   if (i_simu == 2){
  #     feature.fLARS <- c(27, 14, 11, 25, 31, 8, 6, 24, 13)
  #   }else
  #     if (i_simu == 5){
  #       feature.fLARS <- c(11, 13, 21, 22, 25, 30, 32, 1, 5)
  #     }else
  #       if (i_simu == 6){
  #         feature.fLARS <- c(19, 32, 11, 5, 17, 31, 30, 29)
  #       }

    
  # vector of epochs
  epochs_try = rep(256, times = 10)

  numbasis <- 5
  neurons_per_layer <- c(64, 128, 128, 64)

  classcol <- function (col, dim) {
    # converts a category to a column vector of dimension dim
    m <- matrix(0.1, dim, 1)
    m[col, 1] <- 0.9
    m
  }

  class_matrix_to_vector <- function(modeled_output) {
    size <- nrow(modeled_output)
    classvalues <- c(rep(0, size))
    for (i in 1:size) {
      class <- which.max(modeled_output[i, ])
      if (modeled_output[i, class] > 0) {
        classvalues[i] <- class
      } else {
        classvalues[i] = NA
      }
    }
    return(classvalues)
  }

  accuracy <- function(confusion){
    frequent <- 0
    for (i in 1:ncol(confusion)){
      if (colnames(confusion)[i] %in% rownames(confusion)){
        frequent = frequent + confusion[which(rownames(confusion)==colnames(confusion)[i]), i]
      }
    }
    return(frequent / sum(confusion))
  }

  # convert each training.class to a column vector of length nclasses
  # where all entries are zero except for the one at column
  # training.class. Join all the column vectors (for each sample)
  # into a matrix of size training.size by nclasses.
  nclasses <- length(unique(train_y)) # Number of classes
  training.class <- as.numeric(train_y)
  train_class_matrix <- sapply(training.class, classcol, nclasses)
  cdict <- paste0("c",1:nclasses)
  rownames(train_class_matrix) <- cdict
  # Assign variables c1,c2, and etc. to the nrows of the class_matrix.
  
  
  ## FNN with the original features
  # Clearing backend
  K <- backend()
  K$clear_session()
  options(warn=-1)
  # Functional weights & initializations
  func_weights = list()
  fnn_training_plot <- list()

  nvariable <- nsensor
  domain <- list()
  for(tim_len in 1:nvariable) domain[[tim_len]] <- TRange

  # Creating list
  plot_list = list()

  train_fda_full <- train_fda
  test_fda_full <- test_fda
  acc_result <- matrix(nrow = length(epochs_try), ncol = 2)
  for (i in 1:length(epochs_try)) {
    # Getting results
    # i = 1
    fnn_model_full <- fnn.fit(resp = t(train_class_matrix),
                              func_cov = train_fda_full,
                              scalar_cov = NULL,
                              basis_choice = rep(c("fourier"), times = dim(train_fda_full)[3]),
                              num_basis = rep(numbasis, times = dim(train_fda_full)[3]),
                              hidden_layers = length(neurons_per_layer),
                              neurons_per_layer = neurons_per_layer,
                              activations_in_layers = c(rep("relu", times = length(neurons_per_layer))),
                              domain_range = domain,
                              epochs = epochs_try[i],
                              loss_choice = "mse",
                              metric_choice = list("mean_squared_error"),
                              val_split = 0.2,
                              learn_rate = 0.002,
                              patience_param = 50,
                              early_stop = T,
                              print_info = T,
                              batch_size = 32
    )

    # Predicting
    pred_fnn_in_full = fnn.predict(fnn_model_full,
                                   train_fda_full,
                                   scalar_cov = NULL,
                                   basis_choice = rep("fourier", times = dim(train_fda_full)[3]),
                                   num_basis = rep(numbasis, times = dim(train_fda_full)[3]),
                                   domain_range = domain)
    classvalues_in_full <- class_matrix_to_vector(pred_fnn_in_full)

    pred_fnn_out_full = fnn.predict(fnn_model_full,
                                    test_fda_full,
                                    scalar_cov = NULL,
                                    basis_choice = rep("fourier", times = dim(test_fda_full)[3]),
                                    num_basis = rep(numbasis, times = dim(test_fda_full)[3]),
                                    domain_range = domain)
    classvalues_out_full <- class_matrix_to_vector(pred_fnn_out_full)
    ###################
    # Storing Results #
    ###################
    insample_full <- table(train_y, classvalues_in_full)
    training.accuracy <- accuracy(insample_full)
    # cat("training data\n confusion matrix")
    # print(insample_full)
    cat("training accuracy = ", training.accuracy,"\n")

    outsample_full <- table(test_y, classvalues_out_full)
    test.accuracy <- accuracy(outsample_full)
    # cat("test data\n confusion matrix")
    # print(outsample)
    cat("test accuracy = ", test.accuracy,"\n")

    write.csv(insample_full, paste0("confusion_matrix_", simulation[i_simu], "_insample_full_", i, ".csv"))
    write.csv(outsample_full, paste0("confusion_matrix_", simulation[i_simu], "_outsample_full_", i, ".csv"))

    acc_result[i, ] <- c(training.accuracy, test.accuracy)

    # Printing iteration number
    print(paste0("Done Iteration: ", i))


    # Clearning sessions
    K$clear_session()

  }
  colnames(acc_result) <- c("training.accuracy", "test.accuracy")
  write.csv(acc_result, paste0("acc_", simulation[i_simu], "_full.csv"))
  
  
  
  ## FNN with the selected features
  # Clearing backend
  K <- backend()
  K$clear_session()
  options(warn=-1)


  train_fda_selected <- train_fda[,,feature.fb]
  test_fda_selected <- test_fda[,,feature.fb]
  nvariable <- length(feature.fb)
  # train_fda_selected <- train_fda[,,feature.fLARS]
  # test_fda_selected <- test_fda[,,feature.fLARS]
  # nvariable <- length(feature.fLARS)

  ## Functional weights & initializations
  func_weights = list()
  fnn_training_plot <- list()


  domain <- list()
  for(tim_len in 1:nvariable) domain[[tim_len]] <- TRange

  acc_result <- matrix(nrow = length(epochs_try), ncol = 2)
  for (i in 1:length(epochs_try)) {
    # Getting results
    # i = 1
    fnn_model_selected <- fnn.fit(resp = t(train_class_matrix),
                                  func_cov = train_fda_selected,
                                  scalar_cov = NULL,
                                  basis_choice = rep(c("fourier"), times = dim(train_fda_selected)[3]),
                                  num_basis = rep(numbasis, times = dim(train_fda_selected)[3]),
                                  hidden_layers = length(neurons_per_layer),
                                  neurons_per_layer = neurons_per_layer,
                                  activations_in_layers = c(rep("relu", times = length(neurons_per_layer))),
                                  domain_range = domain,
                                  epochs = epochs_try[i],
                                  loss_choice = "mse",
                                  metric_choice = list("mean_squared_error"),
                                  val_split = 0.2,
                                  learn_rate = 0.002,
                                  patience_param = 50,
                                  early_stop = T,
                                  print_info = T,
                                  batch_size = 32
    )

    # Predicting
    pred_fnn_in_selected = fnn.predict(fnn_model_selected,
                                       train_fda_selected,
                                       scalar_cov = NULL,
                                       basis_choice = rep("fourier", times = dim(train_fda_selected)[3]),
                                       num_basis = rep(numbasis, times = dim(train_fda_selected)[3]),
                                       domain_range = domain)
    classvalues_in_selected <- class_matrix_to_vector(pred_fnn_in_selected)

    pred_fnn_out_selected = fnn.predict(fnn_model_selected,
                                        test_fda_selected,
                                        scalar_cov = NULL,
                                        basis_choice = rep("fourier", times = dim(test_fda_selected)[3]),
                                        num_basis = rep(numbasis, times = dim(test_fda_selected)[3]),
                                        domain_range = domain)
    classvalues_out_selected <- class_matrix_to_vector(pred_fnn_out_selected)
    ###################
    # Storing Results #
    ###################
    insample_selected <- table(train_y, classvalues_in_selected)
    training.accuracy <- accuracy(insample_selected)
    # cat("training data\n confusion matrix")
    # print(insample_selected)
    cat("training accuracy = ", training.accuracy,"\n")

    outsample_selected <- table(test_y, classvalues_out_selected)
    test.accuracy <- accuracy(outsample_selected)
    # cat("test data\n confusion matrix")
    # print(outsample_selected)
    cat("test accuracy = ", test.accuracy,"\n")

    write.csv(insample_selected, paste0("confusion_matrix_", simulation[i_simu], "_insample_selected_fofcd_", i, ".csv"))
    write.csv(outsample_selected, paste0("confusion_matrix_", simulation[i_simu], "_outsample_selected_fofcd_", i, ".csv"))

    # write.csv(insample_selected, paste0("confusion_matrix_", simulation[i_simu], "_insample_selected_fLARS_", i, ".csv"))
    # write.csv(outsample_selected, paste0("confusion_matrix_", simulation[i_simu], "_outsample_selected_fLARS_", i, ".csv"))

    acc_result[i, ] <- c(training.accuracy, test.accuracy)

    # Printing iteration number
    print(paste0("Done Iteration: ", i))


    # Clearning sessions
    K$clear_session()

  }
  colnames(acc_result) <- c("training.accuracy", "test.accuracy")
  write.csv(acc_result, paste0("acc_", simulation[i_simu], "_selected_fofcd.csv"))
  # write.csv(acc_result, paste0("acc_", simulation[i_simu], "_selected_fLARS.csv"))
}

