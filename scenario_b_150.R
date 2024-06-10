library(dplyr)
library(Rglpk)
library(tictoc)
library(lpSolveAPI)
library(quadprog)

library("aces")
library("eat")

repl <- 50

simulaciones <- data.frame(
  id = rep(NA, repl),
  nX = rep(NA, repl),
  N = rep(NA, repl),
  degree = rep(NA, repl),
  err_red = rep(NA, repl),
  d = rep(NA, repl),
  metric = rep(NA, repl),
  cluster = rep(NA, repl),
  time = rep(NA, repl),
  mse_aces = rep(NA, repl),
  bias_aces = rep(NA, repl),
  mse_dea = rep(NA, repl),
  bias_dea = rep(NA, repl),
  time_dea = rep(NA, repl)
)

set.seed(314)

for (s in c("B")) {
  for (N in c(150)) {
    for(id in 1:50) {
      
      simulaciones[id, c(1, 2, 3)] <- c(id, 2, N)
      
      while (TRUE) {
        
        # ==
        # Data Generation
        # ==
        
        data <- reffcy (
          DGP = "add_scenario_XnY1",
          parms = list (
            N = N,
            scenario = s
          ))
        
        if (nrow(data) < N) next
        
        # input and output indexes
        if (s %in% c("A", "B")) {
          x <- 1
          y <- 2
        } else if (s %in% c("C", "E")) {
          x <- 1:2
          y <- 3
        } else {
          x <- 1:3
          y <- 4
        }
        
        nX <- length(x)
        nY <- length(y)
        simulaciones[id, 2] <- nX
        
        
        # ==
        # DEA
        # ==
        
        tic()
        
        dea <- aces:::rad_out (
          tech_xmat = as.matrix(data[, x]),
          tech_ymat = as.matrix(data[, y]),
          eval_xmat = as.matrix(data[, x]),
          eval_ymat = as.matrix(data[, y]),
          convexity = TRUE,
          returns = "variable"
        ) * data[, y]
        
        time <- toc()
        
        dif <- dea - data[, "yD"]
        
        # mse | bias | bias absolute
        simulaciones[id, 12] <- round(mean(dif ^ 2), 3)
        simulaciones[id, 13] <- round(mean(dif), 3)
        simulaciones[id, 14] <- unname(time$toc - time$tic)
        
        if (round(sum(dif ^ 2) / N, 3) < 1000) break
      }
      
      print(paste(" ### !!! -- > Réplica número", id, "/", 50, "<-- ¡¡¡ ###"))
      
      # ==
      # ACES
      # ==
      
      grid <- expand.grid (
        err_red = c(0.01),
        d = c(1, 2),
        degree = 1,
        metric = "mae",
        cluster = c(1, 2, 3)
      )
      
      # mse for k-fold cross-validation
      grid$mse1 <- rep(NA, nrow(grid))
      grid$mse2 <- rep(NA, nrow(grid))
      grid$mse3 <- rep(NA, nrow(grid))
      grid$mse4 <- rep(NA, nrow(grid))
      grid$mse5 <- rep(NA, nrow(grid))
      
      # times for k-fold cross-validation
      grid$time1 <- rep(NA, nrow(grid))
      grid$time2 <- rep(NA, nrow(grid))
      grid$time3 <- rep(NA, nrow(grid))
      grid$time4 <- rep(NA, nrow(grid))
      grid$time5 <- rep(NA, nrow(grid))
      
      # mean error in k-fold cross-validation
      grid$mse <- rep(NA, nrow(grid))
      
      # mean time in k-fold cross-validation
      grid$time <- rep(NA, nrow(grid))
      
      # ================ #
      # Cross validation #
      # ================ #
      
      data_shuffle <- data[sample(1:nrow(data)), ]
      
      # number of folds
      kfold <- 5
      
      # size of the fold
      ksize <- floor(N / 5)
      
      # list for each fold
      kindex <- vector("list", kfold)
      
      for (k in 1:kfold) {
        kindex[[k]] <- ((k - 1) * ksize + 1):(k * ksize)
      }
      
      for (k in 1:kfold) {
        
        # training data
        train <- data_shuffle[- kindex[[k]], ]
        
        # test data
        test <- data_shuffle[kindex[[k]] , ]
        
        # clustered data
        for (cl in c(1, 2, 3)) {
          
          if (cl == 1) {
            # assign the unique cluster to all data in training
            train[, paste("cluster_", cl, sep = "")] <- 1
            
            # assign the unique cluster to all data in test
            test[, paste("cluster_", cl, sep = "")] <- 1
            
            next
          }
          
          model_eat <- EAT (
            data = train,
            x = x,
            y = y,
            max.leaves = cl
          )
          
          predict_train_eat <- predict (
            object = model_eat,
            newdata = train,
            x = x
          )
          
          train[, paste("cluster_", cl, sep = "")] <- factor (
            x = predict_train_eat[, 1],
            labels = 1:cl
          )
          
          levels <- sort(unique(predict_train_eat[, 1]))
          
          predict_test_eat <- predict (
            object = model_eat,
            newdata = test,
            x = x
          )
          
          for (p in 1:cl) {
            test[predict_test_eat == levels[p], paste("cluster_", cl, sep = "")] <- p
          }
          
        }
        
        for (k in 1:kfold) {
          
          # training data
          train <- data_shuffle[- kindex[[k]], ]
          
          # test data
          test <- data_shuffle[kindex[[k]], ]
          
          # clustered data
          for (cl in c(1, 2, 3)) {
            
            if (cl == 1) {
              # assign the unique cluster to all data in training
              train[, paste("cluster_", cl, sep = "")] <- 1
              
              # assign the unique cluster to all data in test
              test[, paste("cluster_", cl, sep = "")] <- 1
              
              next
            }
            
            model_eat <- EAT(
              data = train,
              x = x,
              y = y,
              max.leaves = cl
            )
            
            predict_train_eat <- predict(
              object = model_eat,
              newdata = train,
              x = x
            )
            
            train[, paste("cluster_", cl, sep = "")] <- factor(
              x = predict_train_eat[, 1],
              labels = 1:cl
            )
            
            levels <- sort(unique(predict_train_eat[, 1]))
            
            predict_test_eat <- predict(
              object = model_eat,
              newdata = test,
              x = x
            )
            
            for (p in 1:cl) {
              test[predict_test_eat == levels[p], paste("cluster_", cl, sep = "")] <- p
            }
            
          }
          
          for (i in 1:nrow(grid)) {
            
            # number of clusters
            nclusters <- grid[i, "cluster"]
            
            # index of clusters
            cluster_idx <- train[, (nX+nY+1) + nclusters]
            
            # list of mse
            mse_list <- vector("list", nclusters)
            
            for (q in 1:nclusters) {
              mse_list[[q]] <- matrix(0, nrow = 1, ncol = 3)
            }
            
            tic()
            
            for (cl in 1:nclusters) {
              
              # clustered training
              cluster_train <- train[cluster_idx == cl, ]
              
              model <- tryCatch({
                aces(
                  data = cluster_train,
                  x = x,
                  y = y,
                  y_type = "all",
                  model_type = "env",
                  error_type = "add",
                  RF = list(
                    "apply" = FALSE,
                    "sample" = nrow(cluster_train),
                    "models" = 200,
                    "nvars" = 1,
                    "oob_red" = 0.001
                  ),
                  mul_BF = list(
                    "degree" = grid[i, "degree"],
                    "hd_cost" = 0
                  ),
                  metric = grid[i, "metric"],
                  shape = list(
                    "mon" = T,
                    "con" = T,
                    "ori" = F
                  ),
                  nterms = nrow(cluster_train),
                  err_red = grid[i, "err_red"],
                  kn_grid = -1,
                  minspan = -1,
                  endspan = 0,
                  kn_penalty = grid[i, "d"],
                  smoothing = list(
                    "wc" = seq(1, 2, length.out = 5),
                    "wq" = seq(8 / 7, 1.5, length.out = 5)
                  )
                )
              }, error = function(e) {
                message(paste("Error in aces() for cluster", cl, "and grid row", i, ": ", e$message))
                return(NULL)
              })
              
              if (is.null(model)) {
                next
              }
              
              cluster_test <- test[test[, (nX+nY+1) + nclusters] == cl, ]
              
              if (nrow(cluster_test) == 0) {
                next # Saltamos a la siguiente iteración si no hay datos en ese cluster
              }
              
              mse <- matrix(0, nrow = 1, ncol = 3)
              col <- 0
              for (m in c("aces", "aces_cubic", "aces_quintic")) {
                
                col <- col + 1
                
                y_hat_test_cluster <- predict(
                  object = model,
                  newdata = cluster_test,
                  x = x,
                  method = m
                )
                
                dif_test <- y_hat_test_cluster - cluster_test[, c("y")]
                
                mse[1, col] <- round(sum(dif_test[1] ^ 2) / nrow(dif_test), 3)
              }
              
              mse_list[[cl]] <- mse
              
            }
            
            time <- toc()
            
            grid[i, 5 + k] <- sum(
              sapply(
                1:nclusters,
                function(x) (
                  sum(cluster_idx == x) / length(cluster_idx) * (mse_list[[x]][, 1])
                )
              )
            )
            
            grid[i, 10 + k] <- unname(time$toc - time$tic)
          }
        }
        
        
        # mean error in k-fold cross-validation
        grid$mse <- apply(grid[, 06:10], 1, mean)
        
        # mean time in k-fold cross-validation
        grid$time <- apply(grid[, 11:15], 1, mean)
        
        # Best set of hyperparameters
        best_hyp <- grid %>%
          top_n(- 1, mse) %>%
          top_n(- 1, time) %>%
          sample_n(1)
        
        simulaciones[id, 04] <- best_hyp[, "degree"]
        simulaciones[id, 05] <- best_hyp[, "err_red"]
        simulaciones[id, 06] <- best_hyp[, "d"]
        simulaciones[id, 07] <- "mae"
        simulaciones[id, 08] <- best_hyp[, "cluster"]
        
        # EAT with best_cluster
        
        if (best_hyp[, "cluster"] == 1) {
          data[, "cluster"] <- 1
          
          # number of clusters
          nclust <- best_hyp[, "cluster"]
          
          # an aces model for each cluster
          aces_list <- vector("list", nclust)
          
          tic()
          
          for (cl in 1:nclust) {
            
            
            # guardar en la cl-ésima posición un modelo ACES entrenado con el cl-ésimo
            # conjunto de datos
            aces_list[[cl]] <- aces (
              data = data[data[,"cluster"]==cl,],
              x = x,
              y = y,
              y_type = "all",
              model_type = "env",
              error_type = "add",
              RF = list (
                "apply" = FALSE,
                "sample" = nrow(data[data[,"cluster"]==cl,]),
                "models" = 200,
                "nvars" = 1,
                "oob_red" = 0.001
              ),
              mul_BF = list (
                "degree" = simulaciones[id, "degree"],
                "hd_cost" = 0
              ),
              metric = simulaciones[id, "metric"],
              shape = list (
                "mon" = T,
                "con" = T,
                "ori" = F
              ),
              nterms = nrow(data[data[,"cluster"]==cl,]),
              err_red = simulaciones[id, "err_red"],
              kn_grid = - 1,
              minspan = - 1,
              endspan = 0,
              kn_penalty = simulaciones[id, "d"],
              smoothing = list (
                "wc" = seq(1, 2, length.out = 5),
                "wq" = seq(8 / 7, 1.5, length.out = 5)
              )
            )
            
            
          } 
          
          time <- toc()
          
          # time
          simulaciones[id, 09] <- unname(time$toc - time$tic)
          
          # ==
          # aces predictions
          # ==
          
          y_hat_pred <- matrix (
            0, 
            nrow = nrow(data)
          )
          
          for (cl in 1:nclust) {
            
            
            idx_clust <- which(data[,"cluster"]==cl)
            
            for (l in 1:length(idx_clust)){
              
              y_hat_pred[idx_clust[l],] <- predict (
                object = aces_list[[cl]],
                newdata = data[data[,"cluster"]==cl,],
                x = x,
                method = "aces"
              )$y_pred[l]
            }
            
          }
          
          y_hat_pred <- aces:::rad_out (
            tech_xmat = as.matrix(data[, x]),
            tech_ymat = as.matrix(y_hat_pred),
            eval_xmat = as.matrix(data[, x]),
            eval_ymat = as.matrix(data[, y]),
            convexity = TRUE,
            returns = "variable"
          ) * data[, y]
          
          dif <- data[, "yD"] - y_hat_pred
          
          # mse | bias | bias absolute
          simulaciones[id, 10] <- round(mean(dif ^ 2), 3)
          simulaciones[id, 11] <- round(mean(dif), 3)
          #simulaciones[id, cols[j] + 3] <- round(mean(abs(dif)), 3)
          
          next
        }
        
        modelo_eat <- EAT (
          data = data,
          x = x,
          y = y,
          max.leaves = best_hyp[, "cluster"]
        )
        
        # plotEAT(object = modelo_eat)
        
        # predictions
        predictions_EAT <- predict (
          object = modelo_eat,
          newdata = data,
          x = x
        )
        
        predictions_EAT[, 2] <- factor (
          x = predictions_EAT[,1], 
          labels = 1: best_hyp[, "cluster"]
        )
        
        # add number of cluster to data
        data[, "cluster"] <- predictions_EAT[, 2]
        
        # number of clusters
        nclust <- best_hyp[, "cluster"]
        
        # an aces model for each cluster
        aces_list <- vector("list", nclust)
        
        tic()
        
        for (cl in 1:nclust) {
          
          
          # guardar en la cl-ésima posición un modelo ACES entrenado con el cl-ésimo
          # conjunto de datos
          aces_list[[cl]] <- aces (
            data = data[data[,"cluster"]==cl,],
            x = x,
            y = y,
            y_type = "all",
            model_type = "env",
            error_type = "add",
            RF = list (
              "apply" = FALSE,
              "sample" = nrow(data[data[,"cluster"]==cl,]),
              "models" = 200,
              "nvars" = 1,
              "oob_red" = 0.001
            ),
            mul_BF = list (
              "degree" = simulaciones[id, "degree"],
              "hd_cost" = 0
            ),
            metric = simulaciones[id, "metric"],
            shape = list (
              "mon" = T,
              "con" = T,
              "ori" = F
            ),
            nterms = nrow(data[data[,"cluster"]==cl,]),
            err_red = simulaciones[id, "err_red"],
            kn_grid = - 1,
            minspan = - 1,
            endspan = 0,
            kn_penalty = simulaciones[id, "d"],
            smoothing = list (
              "wc" = seq(1, 2, length.out = 5),
              "wq" = seq(8 / 7, 1.5, length.out = 5)
            )
          )
          
          
        } 
        
        time <- toc()
        
        # time
        simulaciones[id, 09] <- unname(time$toc - time$tic)
        
        # ==
        # aces predictions
        # ==
        
        y_hat_pred <- matrix (
          0, 
          nrow = nrow(data)
        )
        
        for (cl in 1:nclust) {
          
          
          idx_clust <- which(data[,"cluster"]==cl)
          
          for (l in 1:length(idx_clust)){
            
            y_hat_pred[idx_clust[l],] <- predict (
              object = aces_list[[cl]],
              newdata = data[data[,"cluster"]==cl,],
              x = x,
              method = "aces"
            )$y_pred[l]
          }
          
        }
        
        y_hat_pred <- aces:::rad_out (
          tech_xmat = as.matrix(data[, x]),
          tech_ymat = as.matrix(y_hat_pred),
          eval_xmat = as.matrix(data[, x]),
          eval_ymat = as.matrix(data[, y]),
          convexity = TRUE,
          returns = "variable"
        ) * data[, y]
        
        dif <- data[, "yD"] - y_hat_pred
        
        # mse | bias | bias absolute
        simulaciones[id, 10] <- round(mean(dif ^ 2), 3)
        simulaciones[id, 11] <- round(mean(dif), 3)
        #simulaciones[id, cols[j] + 3] <- round(mean(abs(dif)), 3)
        
      }
    }
  }
}

aces_mean <- mean(simulaciones[,10])
dea_mean <- mean(simulaciones[,12])

aces_mean
dea_mean


write.csv(simulaciones, "simulaciones_b_nn_x2y1_20_150.csv")
