nsubject <- 27
nsensor <- 10+22
nclass_1 <- 12
nclass_2 <- 17
nclass_B1 <- 8
nclass_B2 <- 9
nclass_3 <- 23
nclass <- 12+17+23
nrep <- 10

## data of class 1
data <- NULL
for (isub in 1:nsubject){
  sensor <- read.csv(paste0('S', isub, '_A1_E1.csv'))
  nsub <- data.frame(nsubject = rep(isub, dim(sensor)[1]))
  data <- rbind(data, data.frame(sensor, nsub))
}

## repetition & stimulus
index <- matrix(rep(0, times = nrep*nclass_1*nsubject), nrow = nclass_1)
for(isub in 1:nsubject){
  for(irep in 1:nrep){
    for(isti in 1:nclass_1){
      index[isti, (isub-1)*nrep+irep] = length(intersect(intersect(which(data$repetition == irep),
                                                                   which(data$stimulus == isti)),
                                                         which(data$nsubject == isub)))
    }
  }
}
ntime <- max(index)

signal <- list()
y <- rep(NA, times = nrep*nclass_1*nsubject)
yclass <- matrix(NA, nrow = nrep*nclass_1*nsubject, ncol = 3)
for(isensor in 1:nsensor){
  signal[[isensor]] <- matrix(NA, nrow = nrep*nclass_1*nsubject, ncol = ntime)
  for(isub in 1:nsubject){
    for(irep in 1:nrep){
      for(isti in 1:nclass_1){
        iindex <- intersect(intersect(which(data$repetition == irep),
                                      which(data$stimulus == isti)),
                            which(data$nsubject == isub))
        signal[[isensor]][(isub-1)*nrep*nclass_1 + (irep-1)*nclass_1 + isti, 1:length(iindex)] <- data[iindex, isensor]
        y[(isub-1)*nrep*nclass_1 + (irep-1)*nclass_1 + isti] <- isti
      }
      yclass[(isub-1)*nrep*nclass_1 + (irep-1)*nclass_1 + 1:nclass_1, 3] <- irep
    }
    yclass[(isub-1)*nrep*nclass_1 + 1:(nclass_1*nrep), 2] <- isub
  }
  write.csv(signal[[isensor]], paste0("Sall_A1_E1_signal_of_sensor_", isensor, ".csv"), row.names = F)
}
yclass[,1] <- y
colnames(yclass) <- c("class", "subject", "repitition")
write.csv(yclass, paste0("Sall_A1_E1_class_of_posture.csv"), row.names = F)


## data of class 2
data <- NULL
for (isub in 1:nsubject){
  sensor <- read.csv(paste0('S', isub, '_A1_E2.csv'))
  nsub <- data.frame(nsubject = rep(isub, dim(sensor)[1]))
  data <- rbind(data, data.frame(sensor, nsub))
}

## repetition & stimulus
index <- matrix(rep(0, times = nrep*nclass_2*nsubject), nrow = nclass_2)
for(isub in 1:nsubject){
  for(irep in 1:nrep){
    for(isti in 1:nclass_2){
      index[isti, (isub-1)*nrep+irep] = length(intersect(intersect(which(data$repetition == irep),
                                                                   which(data$stimulus == isti)),
                                                         which(data$nsubject == isub)))
    }
  }
}
ntime <- max(index)

signal <- list()
y <- rep(NA, times = nrep*nclass_2*nsubject)
yclass <- matrix(NA, nrow = nrep*nclass_2*nsubject, ncol = 3)
for(isensor in 1:nsensor){
  signal[[isensor]] <- matrix(NA, nrow = nrep*nclass_2*nsubject, ncol = ntime)
  for(isub in 1:nsubject){
    for(irep in 1:nrep){
      for(isti in 1:nclass_2){
        iindex <- intersect(intersect(which(data$repetition == irep),
                                      which(data$stimulus == isti)),
                            which(data$nsubject == isub))
        signal[[isensor]][(isub-1)*nrep*nclass_2 + (irep-1)*nclass_2 + isti, 1:length(iindex)] <- data[iindex, isensor]
        y[(isub-1)*nrep*nclass_2 + (irep-1)*nclass_2 + isti] <- isti
      }
      yclass[(isub-1)*nrep*nclass_2 + (irep-1)*nclass_2 + 1:nclass_2, 3] <- irep
    }
    yclass[(isub-1)*nrep*nclass_2 + 1:(nclass_2*nrep), 2] <- isub
  }
  write.csv(signal[[isensor]], paste0("Sall_A1_E2_signal_of_sensor_", isensor, ".csv"), row.names = F)
}
yclass[,1] <- y
colnames(yclass) <- c("class", "subject", "repitition")
write.csv(yclass, paste0("Sall_A1_E2_class_of_posture.csv"), row.names = F)


## data of class 3
data <- NULL
for (isub in 1:nsubject){
  sensor <- read.csv(paste0('S', isub, '_A1_E3.csv'))
  nsub <- data.frame(nsubject = rep(isub, dim(sensor)[1]))
  data <- rbind(data, data.frame(sensor, nsub))
}

## repetition & stimulus
index <- matrix(rep(0, times = nrep*nclass_3*nsubject), nrow = nclass_3)
for(isub in 1:nsubject){
  for(irep in 1:nrep){
    for(isti in 1:nclass_3){
      index[isti, (isub-1)*nrep+irep] = length(intersect(intersect(which(data$repetition == irep),
                                                                   which(data$stimulus == isti)),
                                                         which(data$nsubject == isub)))
    }
  }
}
ntime <- max(index)

signal <- list()
y <- rep(NA, times = nrep*nclass_3*nsubject)
yclass <- matrix(NA, nrow = nrep*nclass_3*nsubject, ncol = 3)
for(isensor in 1:nsensor){
  signal[[isensor]] <- matrix(NA, nrow = nrep*nclass_3*nsubject, ncol = ntime)
  for(isub in 1:nsubject){
    for(irep in 1:nrep){
      for(isti in 1:nclass_3){
        iindex <- intersect(intersect(which(data$repetition == irep),
                                      which(data$stimulus == isti)),
                            which(data$nsubject == isub))
        signal[[isensor]][(isub-1)*nrep*nclass_3 + (irep-1)*nclass_3 + isti, 1:length(iindex)] <- data[iindex, isensor]
        y[(isub-1)*nrep*nclass_3 + (irep-1)*nclass_3 + isti] <- isti
      }
      yclass[(isub-1)*nrep*nclass_3 + (irep-1)*nclass_3 + 1:nclass_3, 3] <- irep
    }
    yclass[(isub-1)*nrep*nclass_3 + 1:(nclass_3*nrep), 2] <- isub
  }
  write.csv(signal[[isensor]], paste0("Sall_A1_E3_signal_of_sensor_", isensor, ".csv"), row.names = F)
}
yclass[,1] <- y
colnames(yclass) <- c("class", "subject", "repitition")
write.csv(yclass, paste0("Sall_A1_E3_class_of_posture.csv"), row.names = F)


## data of all classes
data <- NULL
for (isub in 1:nsubject){
  sensor <- read.csv(paste0('S', isub, '_A1_E1.csv'))
  nsub <- data.frame(nsubject = rep(isub, dim(sensor)[1]))
  data <- rbind(data, data.frame(sensor, nsub))
  
  sensor <- read.csv(paste0('S', isub, '_A1_E2.csv'))
  nsub <- data.frame(nsubject = rep(isub, dim(sensor)[1]))
  stimulus_n0 <- which(sensor$stimulus != 0)
  sensor$stimulus[stimulus_n0] <- sensor$stimulus[stimulus_n0] + nclass_1
  restimulus_n0 <- which(sensor$restimulus != 0)
  sensor$restimulus[restimulus_n0] <- sensor$restimulus[restimulus_n0] + nclass_1
  data <- rbind(data, data.frame(sensor, nsub))
  
  sensor <- read.csv(paste0('S', isub, '_A1_E3.csv'))
  nsub <- data.frame(nsubject = rep(isub, dim(sensor)[1]))
  stimulus_n0 <- which(sensor$stimulus != 0)
  sensor$stimulus[stimulus_n0] <- sensor$stimulus[stimulus_n0] + nclass_1 + nclass_2
  restimulus_n0 <- which(sensor$restimulus != 0)
  sensor$restimulus[restimulus_n0] <- sensor$restimulus[restimulus_n0] + nclass_1 + nclass_2
  data <- rbind(data, data.frame(sensor, nsub))
}

# index = intersect(intersect(which(data$repetition == 1), which(data$stimulus == 1)),
#                   which(data$nsubject == 1))
# plot(data[index, 1], type = "l", col = "grey")
# lines(data[index, 2], type = "l", col = "red")
# lines(data[index, 3], type = "l", col = "blue")
# lines(data[index, 4], type = "l", col = "green")
# lines(data[index, 5], type = "l", col = "yellow")


## repetition & stimulus
index <- matrix(rep(0, times = nrep*nclass*nsubject), nrow = nclass)
for(isub in 1:nsubject){
  for(irep in 1:nrep){
    for(isti in 1:nclass){
      index[isti, (isub-1)*nrep+irep] = length(intersect(intersect(which(data$repetition == irep),
                                                     which(data$stimulus == isti)),
                                           which(data$nsubject == isub)))
    }
  }
}
ntime <- max(index)

signal <- list()
y <- rep(NA, times = nrep*nclass*nsubject)
yclass <- matrix(NA, nrow = nrep*nclass*nsubject, ncol = 3)
for(isensor in 1:nsensor){
  signal[[isensor]] <- matrix(NA, nrow = nrep*nclass*nsubject, ncol = ntime)
  for(isub in 1:nsubject){
    for(irep in 1:nrep){
      for(isti in 1:nclass){
        iindex <- intersect(intersect(which(data$repetition == irep),
                                      which(data$stimulus == isti)),
                            which(data$nsubject == isub))
        signal[[isensor]][(isub-1)*nrep*nclass + (irep-1)*nclass + isti, 1:length(iindex)] <- data[iindex, isensor]
        y[(isub-1)*nrep*nclass + (irep-1)*nclass + isti] <- isti
        # print((irep-1)*nclass + isti)
        # print((isub-1)*nrep*nclass + (irep-1)*nclass + isti)
      }
      yclass[(isub-1)*nrep*nclass + (irep-1)*nclass + 1:nclass, 3] <- irep
    }
    yclass[(isub-1)*nrep*nclass + 1:(nclass*nrep), 2] <- isub
  }
  write.csv(signal[[isensor]], paste0("Sall_A1_Eall_signal_of_sensor_", isensor, ".csv"), row.names = F)
}
yclass[,1] <- y
colnames(yclass) <- c("class", "subject", "repitition")
write.csv(yclass, paste0("Sall_A1_Eall_class_of_posture.csv"), row.names = F)


