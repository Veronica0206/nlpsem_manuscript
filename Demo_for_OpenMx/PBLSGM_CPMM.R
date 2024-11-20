getCPMM_PBLSGM <- function(dat, T_records, nClass, traj_var, t_var, clus_cov, res_ratio = list(rep(4, 2), rep(4, 2)), rho = list(rep(0.1, 2), rep(0.1, 2)),
                           btw_rho = list(rep(0.1, 1), rep(0.1, 1)), btw_res = list(rep(0.3, 1), rep(0.3, 1)), starts = NA, prop_starts = rep(0.5, 2),
                           loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA){
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  dat$ONE <- 1
  if (sum(prop_starts) != 1){
    stop("The sum of all proportion components should be 1!")
  }
  else{
    nT <- rep(0, length(traj_var))
    for (traj in 1:length(traj_var)){
      nT[traj] <- length(T_records[[traj]])
    }
    if (any(is.na(starts))){
      starts <- list()
      FMM <- getFMM_PBLSGM(dat = dat, T_records = T_records, nClass = nClass, traj_var = traj_var, t_var = t_var, res_ratio = res_ratio,
                           rho = rho, btw_rho = btw_rho, btw_res = btw_res, prop_starts = prop_starts, original = F)
      model.para <- summary(FMM)$parameters[, c(1, 5, 6)]
      for (k in 1:nClass){
        uni_mean0.s <- uni_var0.s <- uni_residual <- list()
        mug <- model.para$Estimate[c(grep(paste0("c", k, "muknot"), model.para$name))]
        for (traj in 1:length(traj_var)){
          uni_mean0.s[[length(uni_mean0.s) + 1]] <- c(mxEvalByName(paste0("c", k, "mean_s", traj_var[traj]), model = FMM@submodels[[k]]), mug[traj])
          uni_var0.s[[length(uni_var0.s) + 1]] <- mxEvalByName(paste0("c", k, "psi_s", traj_var[traj]), model = FMM@submodels[[k]])
          uni_residual[[length(uni_residual) + 1]] <- FMM@submodels[[k]]$S$values[((traj - 1) * nT[traj] + 1), ((traj - 1) * nT[traj] + 1)]
        }
        multi_var0.s <- as.matrix(Matrix::bdiag(uni_var0.s))
        multi_residuals <- as.matrix(Matrix::bdiag(uni_residual))
        for (traj_i in 1:(length(traj_var) - 1)){
          for (traj_j in traj_i:(length(traj_var) - 1)){
            multi_var0.s[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
              mxEvalByName(paste0("c", k, "psi_s", traj_var[traj_i], traj_var[traj_j + 1]), model = FMM@submodels[[k]])
            multi_residuals[traj_i, (traj_j + 1)] <- FMM@submodels[[k]]$S$values[((traj_i - 1) * nT[traj_i] + 1), (traj_j * nT[traj_j] + 1)]
          }
        }
        starts[[length(starts) + 1]] <- list(uni_mean0.s, multi_var0.s, multi_residuals)
      }
    }
    FMM_post <- getPosterior(model = FMM, classProbs = FMM$weights$values, round = 4)
    dat$label <- apply(FMM_post, 1, which.max)
    dat_nnet <- dat[, c("label", clus_cov)]
    mod_nnet <- nnet::multinom(label ~ ., data = dat_nnet)
    beta_starts <- (as.matrix(rbind(rep(0, length(clus_cov) + 1), summary(mod_nnet)$coefficient)))
    rownames(beta_starts) <- colnames(beta_starts) <- NULL
  }
  ### Define manifest variables
  traj_list <- list()
  for (traj in 1:length(traj_var)){
    traj_list[[length(traj_list) + 1]] <- paste0(traj_var[traj], T_records[[traj]])
  }
  manifests <- unlist(traj_list)
  ### Define latent variables
  latents <- paste0(rep(c("eta0s", "eta1s", "eta2s"), length(traj_var)), rep(traj_var, each = 3))
  class.list <- list()
  for(k in 1:nClass){
    outDef <- outLoads1 <- outLoads2 <- outDef_L <- outLoads1_L <- outLoads2_L <- list()
    loadings1 <- loadings2 <- loadings3 <- gf_mean <- gamma_mean <- residual_var <- residual_cor <- list()
    func_L <- grad_L <- mean_s_L <- psi_s_L <- psi_btw_s_L <- mean_L <- psi_L <- psi_btw_L <- list()
    gf_var_label <- gf_cov_label <- list()
    for (traj in 1:length(traj_var)){
      for (j in T_records[[traj]]){
        outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[traj], j),
                                name = paste0(traj_var[traj], "t", j))
        outLoads1[[j]] <- mxAlgebraFromString(paste0(traj_var[traj], "t", j, " - c", k, "mug", traj_var[traj]),
                                              name = paste0("c", k, "L1", j, traj_var[traj]))
        outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(", traj_var[traj], "t", j, " - c", k, "mug", traj_var[traj], ")"),
                                              name = paste0("c", k, "L2", j, traj_var[traj]))
      }
      outDef_L[[length(outDef_L) + 1]] <- outDef
      outLoads1_L[[length(outLoads1_L) + 1]] <- outLoads1
      outLoads2_L[[length(outLoads2_L) + 1]] <- outLoads2
      outDef <- outLoads1 <- outLoads2 <- list()
      #### Define factor loadings from latent variables to manifests
      loadings1[[length(loadings1) + 1]] <- mxPath(from = paste0("eta0s", traj_var[traj]), to = traj_list[[traj]],
                                                   arrows = 1, free = F, values = 1)
      loadings2[[length(loadings2) + 1]]  <- mxPath(from = paste0("eta1s", traj_var[traj]), to = traj_list[[traj]],
                                                    arrows = 1, free = F, values = 0,
                                                    labels = paste0("c", k, "L1", T_records[[traj]], traj_var[traj], "[1,1]"))
      loadings3[[length(loadings3) + 1]]  <- mxPath(from = paste0("eta2s", traj_var[traj]), to = traj_list[[traj]],
                                                    arrows = 1, free = F, values = 0,
                                                    labels = paste0("c", k, "L2", T_records[[traj]], traj_var[traj], "[1,1]"))
      #### Define means of outcome-specific growth factors
      gf_mean[[length(gf_mean) + 1]] <- mxPath(from = "one", to = paste0(c("eta0s", "eta1s", "eta2s"), traj_var[traj]),
                                               arrows = 1, free = T, values = starts[[k]][[1]][[traj]][1:3],
                                               labels = paste0("c", k, c("mueta0s", "mueta1s", "mueta2s"), traj_var[traj]))
      #### Add additional parameter and constraints
      gamma_mean[[length(gamma_mean) + 1]] <- mxMatrix("Full", 1, 1, free = T, values = starts[[k]][[1]][[traj]][4],
                                                       labels = paste0("c", k, "muknot_", traj_var[traj]),
                                                       name = paste0("c", k, "mug", traj_var[traj]))
      #### Define the variances of residuals
      residual_var[[length(residual_var) + 1]] <- mxPath(from = traj_list[[traj]], to = traj_list[[traj]], arrows = 2, free = T,
                                                         values = starts[[k]][[3]][traj, traj],
                                                         labels = paste0("c", k, "residuals", traj_var[traj]))
      #### Define transformed function
      func_L[[length(func_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "mug", traj_var[traj], ", c", k, "mug", traj_var[traj], "), ",
                                                                 "cbind(0, 1, -1), ", "cbind(0, 1, 1))"),
                                                          name = paste0("c", k, "func", traj_var[traj]))
      grad_L[[length(grad_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -c", k, "mug", traj_var[traj], ", c", k, "mug", traj_var[traj], "), ",
                                                                 "cbind(0, 1, -1), ", "cbind(0, 1, 1))"),
                                                          name = paste0("c", k, "grad", traj_var[traj]))
      #### Define the outcome-specific growth factor mean vector in the reparameterized framework
      mean_s_L[[length(mean_s_L) + 1]] <- mxAlgebraFromString(paste0("rbind(c", k, "mueta0s", traj_var[traj],
                                                                     ", c", k, "mueta1s", traj_var[traj],
                                                                     ", c", k, "mueta2s", traj_var[traj], ")"),
                                                              name = paste0("c", k, "mean_s", traj_var[traj]))
      #### Define the outcome-specific growth factor var-cov matrix in the reparameterized framework
      psi_s_L[[length(psi_s_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "psi0s0s", traj_var[traj], traj_var[traj],
                                                                   ", c", k, "psi0s1s", traj_var[traj], traj_var[traj],
                                                                   ", c", k, "psi0s2s", traj_var[traj], traj_var[traj], "), ",
                                                                   "cbind(c", k, "psi0s1s", traj_var[traj], traj_var[traj],
                                                                   ", c", k, "psi1s1s", traj_var[traj], traj_var[traj],
                                                                   ", c", k, "psi1s2s", traj_var[traj], traj_var[traj], "), ",
                                                                   "cbind(c", k, "psi0s2s", traj_var[traj], traj_var[traj],
                                                                   ", c", k, "psi1s2s", traj_var[traj], traj_var[traj],
                                                                   ", c", k, "psi2s2s", traj_var[traj], traj_var[traj], "))"),
                                                            name = paste0("c", k, "psi_s", traj_var[traj]))
      #### Define the outcome-specific growth factor mean vector in the original framework
      mean_L[[length(mean_L) + 1]] <- mxAlgebraFromString(paste0("c", k, "func", traj_var[traj], " %*% c", k, "mean_s", traj_var[traj]),
                                                          name = paste0("c", k, "mean", traj_var[traj]))
      #### Define the outcome-specific growth factor var-cov matrix in the original framework
      psi_L[[length(psi_L) + 1]] <- mxAlgebraFromString(paste0("c", k, "grad", traj_var[traj], " %*% c", k, "psi_s", traj_var[traj],
                                                               " %*% t(c", k, "grad", traj_var[traj], ")"),
                                                        name = paste0("c", k, "psi", traj_var[traj]))
      #### Define var-cov of outcome-specific growth factors
      gf_var_label[[length(gf_var_label) + 1]] <- matrix(paste0("c", k, "psi", c("0s0s", "0s1s", "0s2s", "1s0s", "1s1s", "1s2s", "2s0s", "2s1s", "2s2s"),
                                                                traj_var[traj], traj_var[traj]), byrow = T, nrow = 3, ncol = 3)
    }
    for (traj_i in 1:(length(traj_var) - 1)){
      for (traj_j in traj_i:(length(traj_var) - 1)){
        #### Define the between-outcome growth factor var-cov matrix in the reparameterized framework
        psi_btw_s_L[[length(psi_btw_s_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(c", k, "psi0s0s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                             ", c", k, "psi0s1s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                             ", c", k, "psi0s2s", traj_var[traj_i], traj_var[traj_j + 1], "),",
                                                                             "cbind(c", k, "psi1s0s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                             ", c", k, "psi1s1s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                             ", c", k, "psi1s2s", traj_var[traj_i], traj_var[traj_j + 1], "), ",
                                                                             "cbind(c", k, "psi2s0s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                             ", c", k, "psi2s1s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                             ", c", k, "psi2s2s", traj_var[traj_i], traj_var[traj_j + 1],"))"),
                                                                      name = paste0("c", k, "psi_s", traj_var[traj_i], traj_var[traj_j + 1]))
        #### Define the between-outcome growth factor var-cov matrix in the reparameterized framework
        psi_btw_L[[length(psi_btw_L) + 1]] <- mxAlgebraFromString(paste0("c", k, "grad", traj_var[traj_i], " %*% c", k, "psi_s",
                                                                         traj_var[traj_i], traj_var[traj_j + 1],
                                                                         " %*% t(c", k, "grad", traj_var[traj_i + 1], ")"),
                                                                  name = paste0("c", k, "psi", traj_var[traj_i], traj_var[traj_j + 1]))
        #### Define the covariances of residuals
        if (setequal(substr(traj_list[[traj_i]], 2, 2), substr(traj_list[[traj_j + 1]], 2, 2))){
          residual_cor[[length(residual_cor) + 1]] <- mxPath(from = traj_list[[traj_i]], to = traj_list[[traj_j + 1]],
                                                             arrows = 2, free = T, values = starts[[k]][[3]][traj_i, traj_j + 1],
                                                             labels = paste0("c", k, "residuals", traj_var[traj_i], traj_var[traj_j + 1]))
        }
        else{
          T_common <- Reduce(intersect, list(substr(traj_list[[traj_i]], 2, 2), substr(traj_list[[traj_j + 1]], 2, 2)))
          residual_cor[[length(residual_cor) + 1]] <- mxPath(from = paste0(traj_var[traj_i], T_common),
                                                             to = paste0(traj_var[traj_j + 1], T_common),
                                                             arrows = 2, free = T, values = starts[[k]][[3]][traj_i, traj_j + 1],
                                                             labels = paste0("c", k, "residuals", traj_var[traj_i], traj_var[traj_j + 1]))
        }
        gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0("c", k, "psi", c("0s0s", "0s1s", "0s2s", "1s0s", "1s1s", "1s2s", "2s0s", "2s1s", "2s2s"),
                                                             traj_var[traj_i], traj_var[traj_j + 1]), byrow = T, nrow = 3, ncol = 3)
      }
    }
    multi_label <- matrix(NA, nrow = length(latents), ncol = length(latents))
    for (traj in 1:length(traj_var)){
      multi_label[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <- gf_var_label[[traj]]
    }
    for (traj_i in 1:(length(traj_var) - 1)){
      for (traj_j in traj_i:(length(traj_var) - 1)){
        multi_label[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <- gf_cov_label[[traj_i + traj_j - 1]]
      }
    }
    ### Create a mxModel object
    class.list[[length(class.list) + 1]] <- mxModel(paste0("Class", k), type = "RAM",
                                                    manifestVars = manifests, latentVars = latents,
                                                    mxData(observed = dat, type = "raw"),
                                                    #### Define var-cov matrix of latent variables
                                                    mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                                                           free = T, values = t(starts[[k]][[2]])[row(t(starts[[k]][[2]])) >= col(t(starts[[k]][[2]]))],
                                                           labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                                                    outDef_L, outLoads1_L, outLoads2_L, loadings1, loadings2, loadings3,
                                                    gf_mean, gamma_mean, residual_var, residual_cor, func_L, grad_L,
                                                    mean_s_L, psi_s_L, psi_btw_s_L, mean_L, psi_L, psi_btw_L,
                                                    mxFitFunctionML(vector = T))
  }
  ### Make the class proportion matrix, fixing one parameter at a non-zero constant (one)
  classBeta <- mxMatrix(type = "Full", nrow = nClass, ncol = length(clus_cov) + 1,
                        free = rep(c(F, rep(T, nClass - 1)), length(clus_cov) + 1), values = beta_starts,
                        labels = paste0("beta", rep(1:nClass), rep(0:length(clus_cov), each = nClass)),
                        name = "classbeta")
  classPV <- mxMatrix(nrow = length(clus_cov) + 1, ncol = 1, labels = paste0("data.", c("ONE", clus_cov)), values = 1, name = "weightsV")
  classP <- mxAlgebra(classbeta %*% weightsV, name = "weights")
  algebraObjective <- mxExpectationMixture(paste0("Class", 1:nClass),
                                           weights = "weights", scale = "softmax")
  objective <- mxFitFunctionML()
  model_mx <- mxModel("Cluster Predictor Mixture Model, Parallel Bilinear Spline Growth Model with Fixed Knots",
                      mxData(observed = dat, type = "raw"), class.list, classBeta,
                      classPV, classP, algebraObjective, objective)
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if (original){
    model.para <- summary(model)$parameters[, c(1, 5, 6)]
    model.est <- model.se <- est <- list()
    for (k in 1:nClass){
      mean_est <- mean_se <- psi_est <- psi_se <- psi_btw_est <- psi_btw_se <- outcome_est <- outcome_se <- btw_est <- btw_se <- list()
      mug <- model.para$Estimate[c(grep(paste0("c", k, "muknot"), model.para$name))]
      mug_se <- model.para$Std.Error[c(grep(paste0("c", k, "muknot"), model.para$name))]
      for (traj in 1:length(traj_var)){
        mean_est[[length(mean_est) + 1]] <- c(mxEvalByName(paste0("c", k, "mean", traj_var[traj]), model = model@submodels[[k]]), mug[traj])
        mean_se[[length(mean_se) + 1]] <- c(mxSE(paste0("Class", k, ".c", k, "mean", traj_var[traj]), model, forceName = T), mug_se[traj])
        psi_est[[length(psi_est) + 1]] <- mxEvalByName(paste0("c", k, "psi", traj_var[traj]), model = model@submodels[[k]])
        psi_se[[length(psi_se) + 1]] <- mxSE(paste0("Class", k, ".c", k, "psi", traj_var[traj]), model, forceName = T)
        outcome_est[[length(outcome_est) + 1]] <- c(unlist(mean_est[[traj]]), psi_est[[traj]][row(psi_est[[traj]]) >= col(psi_est[[traj]])])
        outcome_se[[length(outcome_se) + 1]] <- c(unlist(mean_se[[traj]]), psi_se[[traj]][row(psi_se[[traj]]) >= col(psi_se[[traj]])])
      }
      for (traj in 1:(length(traj_var) - 1)){
        psi_btw_est[[length(psi_btw_est) + 1]] <- mxEvalByName(paste0("c", k, "psi", traj_var[traj], traj_var[traj + 1]), model = model@submodels[[k]])
        psi_btw_se[[length(psi_btw_se) + 1]] <- mxSE(paste0("Class", k, ".c", k, "psi", traj_var[traj], traj_var[traj + 1]), model, forceName = T)
        btw_est[[length(btw_est) + 1]] <- unlist(c(psi_btw_est))
        btw_se[[length(btw_se) + 1]] <- unlist(c(psi_btw_se))
      }
      model.est[[length(model.est) + 1]] <- round(c(unlist(outcome_est), unlist(btw_est), model.para[grep(paste0("c", k, "residual"), model.para$name), 2]), 4)
      model.se[[length(model.se) + 1]] <- round(c(unlist(outcome_se), unlist(btw_se), model.para[grep(paste0("c", k, "residual"), model.para$name), 3]), 4)
      est[[length(est) + 1]] <- data.frame(Name = paste0("c", k, paraNames),
                                           Estimate = model.est[[k]], SE = model.se[[k]])
    }
    est.beta <- data.frame(Name = paste0("beta", rep(2:nClass, length(clus_cov)), rep(0:2, each = nClass - 1)),
                           Estimate = c(mxEval(classbeta, model)[-1, ]),
                           SE = c(mxSE(classbeta, model)[-1, ]))
    estimate_out <- rbind(do.call(rbind.data.frame, est), est.beta)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}
