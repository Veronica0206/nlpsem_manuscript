getPBLSGM_Fixed <- function(dat, T_records, traj_var, t_var, res_ratio = rep(4, 2), rho = rep(0.1, 2), btw_rho = rep(0.1, 1), btw_res = rep(0.1, 1),
                            starts = NA, loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA){
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  if (any(is.na(starts))){
    uni_residual <- list()
    uni_mean0.s <- uni_var0 <- uni_cov0 <- uni_grad0 <- uni_var0.s <- uni_cov0.s <- list()
    for (traj in 1:length(traj_var)){
      BLSGM_F <- getBLSGM_Fixed(dat = dat, T_records = T_records[[traj]], traj_var = traj_var[traj], t_var = t_var[traj],
                                res_ratio = res_ratio[traj], rho = rho[traj], original = F)
      ### Mean vector of growth factors (in the reparameterized framework)
      uni_mean0.s[[length(uni_mean0.s) + 1]] <- c(BLSGM_F$mean_s$result, BLSGM_F$mug$values)
      ### Residual variance
      uni_residual[[length(uni_residual) + 1]] <- BLSGM_F$output$estimate[1]
      ### var-cov matrix (in the original framework)
      ###################
      uni_var0[[length(uni_var0) + 1]] <- BLSGM_F$psi$result
      uni_grad0[[length(uni_grad0) + 1]] <- matrix(c(1, BLSGM_F$mug$values, 0,
                                                     0, 0.5, 0.5,
                                                     0, -0.5, 0.5), nrow = 3, byrow = T)
    }
    multi_var0 <- as.matrix(Matrix::bdiag(uni_var0))
    multi_grad0 <- as.matrix(Matrix::bdiag(uni_grad0))
    multi_residuals <- as.matrix(Matrix::bdiag(uni_residual))
    for (traj_i in 1:(length(traj_var) - 1)){
      for (traj_j in traj_i:(length(traj_var) - 1)){
        multi_var0[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
          btw_rho[traj_i] * sqrt(ifelse(diag(uni_var0[[traj_i]]) %*% t(diag(uni_var0[[traj_j + 1]])) >= 0,
                                        diag(uni_var0[[traj_i]]) %*% t(diag(uni_var0[[traj_j + 1]])), NA))
        multi_residuals[traj_i, (traj_j + 1)] <- btw_res[traj_i] * sqrt(uni_residual[[traj_i]] * uni_residual[[traj_j + 1]])
      }
    }
    multi_var0.s <- multi_grad0 %*% multi_var0 %*% t(multi_grad0)
    starts <- list(uni_mean0.s, multi_var0.s, multi_residuals)
  }
  ### Define manifest variables
  traj_list <- list()
  for (traj in 1:length(traj_var)){
    traj_list[[length(traj_list) + 1]] <- paste0(traj_var[traj], T_records[[traj]])
  }
  manifests <- unlist(traj_list)
  ### Define latent variables
  latents <- paste0(rep(c("eta0s", "eta1s", "eta2s"), length(traj_var)), rep(traj_var, each = 3))
  outDef <- outLoads1 <- outLoads2 <- outDef_L <- outLoads1_L <- outLoads2_L <- list()
  loadings1 <- loadings2 <- loadings3 <- gf_mean <- gamma_mean <- residual_var <- residual_cor <- list()
  func_L <- grad_L <- mean_s_L <- psi_s_L <- psi_btw_s_L <- mean_L <- psi_L <- psi_btw_L <- list()
  gf_var_label <- gf_cov_label <- list()
  for (traj in 1:length(traj_var)){
    for (j in T_records[[traj]]){
      outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[traj], j),
                              name = paste0(traj_var[traj], "t", j))
      outLoads1[[j]] <- mxAlgebraFromString(paste0(traj_var[traj], "t", j, " - mug", traj_var[traj]),
                                            name = paste0("L1", j, traj_var[traj]))
      outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(", traj_var[traj], "t", j, " - mug", traj_var[traj], ")"),
                                            name = paste0("L2", j, traj_var[traj]))
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
                                                  labels = paste0("L1", T_records[[traj]], traj_var[traj], "[1,1]"))
    loadings3[[length(loadings3) + 1]]  <- mxPath(from = paste0("eta2s", traj_var[traj]), to = traj_list[[traj]],
                                                  arrows = 1, free = F, values = 0,
                                                  labels = paste0("L2", T_records[[traj]], traj_var[traj], "[1,1]"))
    #### Define means of outcome-specific growth factors
    gf_mean[[length(gf_mean) + 1]] <- mxPath(from = "one", to = paste0(c("eta0s", "eta1s", "eta2s"), traj_var[traj]),
                                             arrows = 1, free = T, values = starts[[1]][[traj]][1:3],
                                             labels = paste0(c("mueta0s", "mueta1s", "mueta2s"), traj_var[traj]))
    #### Add additional parameter and constraints
    gamma_mean[[length(gamma_mean) + 1]] <- mxMatrix("Full", 1, 1, free = T, values = starts[[1]][[traj]][4],
                                                     labels = paste0("muknot_", traj_var[traj]),
                                                     name = paste0("mug", traj_var[traj]))
    #### Define the variances of residuals
    residual_var[[length(residual_var) + 1]] <- mxPath(from = traj_list[[traj]], to = traj_list[[traj]], arrows = 2, free = T,
                                                       values = starts[[3]][traj, traj],
                                                       labels = paste0("residuals", traj_var[traj]))
    #### Define transformed function
    func_L[[length(func_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -mug", traj_var[traj], ", mug", traj_var[traj], "), ",
                                                               "cbind(0, 1, -1), ", "cbind(0, 1, 1))"),
                                                        name = paste0("func", traj_var[traj]))
    grad_L[[length(grad_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(1, -mug", traj_var[traj], ", mug", traj_var[traj], "), ",
                                                               "cbind(0, 1, -1), ", "cbind(0, 1, 1))"),
                                                        name = paste0("grad", traj_var[traj]))
    #### Define the outcome-specific growth factor mean vector in the reparameterized framework
    mean_s_L[[length(mean_s_L) + 1]] <- mxAlgebraFromString(paste0("rbind(mueta0s", traj_var[traj], ", mueta1s", traj_var[traj],
                                                                   ", mueta2s", traj_var[traj], ")"),
                                                            name = paste0("mean_s", traj_var[traj]))
    #### Define the outcome-specific growth factor var-cov matrix in the reparameterized framework
    psi_s_L[[length(psi_s_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(psi0s0s", traj_var[traj], traj_var[traj],
                                                                 ", psi0s1s", traj_var[traj], traj_var[traj],
                                                                 ", psi0s2s", traj_var[traj], traj_var[traj], "), ",
                                                                 "cbind(psi0s1s", traj_var[traj], traj_var[traj],
                                                                 ", psi1s1s", traj_var[traj], traj_var[traj],
                                                                 ", psi1s2s", traj_var[traj], traj_var[traj], "), ",
                                                                 "cbind(psi0s2s", traj_var[traj], traj_var[traj],
                                                                 ", psi1s2s", traj_var[traj], traj_var[traj],
                                                                 ", psi2s2s", traj_var[traj], traj_var[traj], "))"),
                                                          name = paste0("psi_s", traj_var[traj]))
    #### Define the outcome-specific growth factor mean vector in the original framework
    mean_L[[length(mean_L) + 1]] <- mxAlgebraFromString(paste0("func", traj_var[traj], " %*% mean_s", traj_var[traj]),
                                                        name = paste0("mean", traj_var[traj]))
    #### Define the outcome-specific growth factor var-cov matrix in the original framework
    psi_L[[length(psi_L) + 1]] <- mxAlgebraFromString(paste0("grad", traj_var[traj], " %*% psi_s", traj_var[traj],
                                                             " %*% t(grad", traj_var[traj], ")"),
                                                      name = paste0("psi", traj_var[traj]))
    #### Define var-cov of outcome-specific growth factors
    gf_var_label[[length(gf_var_label) + 1]] <- matrix(paste0("psi", c("0s0s", "0s1s", "0s2s", "1s0s", "1s1s", "1s2s", "2s0s", "2s1s", "2s2s"),
                                                              traj_var[traj], traj_var[traj]), byrow = T, nrow = 3, ncol = 3)
  }
  for (traj_i in 1:(length(traj_var) - 1)){
    for (traj_j in traj_i:(length(traj_var) - 1)){
      #### Define the between-outcome growth factor var-cov matrix in the reparameterized framework
      psi_btw_s_L[[length(psi_btw_s_L) + 1]] <- mxAlgebraFromString(paste0("rbind(cbind(psi0s0s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                           ", psi0s1s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                           ", psi0s2s", traj_var[traj_i], traj_var[traj_j + 1], "),",
                                                                           "cbind(psi1s0s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                           ", psi1s1s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                           ", psi1s2s", traj_var[traj_i], traj_var[traj_j + 1], "), ",
                                                                           "cbind(psi2s0s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                           ", psi2s1s", traj_var[traj_i], traj_var[traj_j + 1],
                                                                           ", psi2s2s", traj_var[traj_i], traj_var[traj_j + 1],"))"),
                                                                    name = paste0("psi_s", traj_var[traj_i], traj_var[traj_j + 1]))
      #### Define the between-outcome growth factor var-cov matrix in the reparameterized framework
      psi_btw_L[[length(psi_btw_L) + 1]] <- mxAlgebraFromString(paste0("grad", traj_var[traj_i], " %*% psi_s",
                                                                       traj_var[traj_i], traj_var[traj_j + 1],
                                                                       " %*% t(grad", traj_var[traj_i + 1], ")"),
                                                                name = paste0("psi", traj_var[traj_i], traj_var[traj_j + 1]))
      #### Define the covariances of residuals
      if (setequal(substr(traj_list[[traj_i]], 2, 2), substr(traj_list[[traj_j + 1]], 2, 2))){
        residual_cor[[length(residual_cor) + 1]] <- mxPath(from = traj_list[[traj_i]], to = traj_list[[traj_j + 1]],
                                                           arrows = 2, free = T, values = starts[[3]][traj_i, traj_j + 1],
                                                           labels = paste0("residuals", traj_var[traj_i], traj_var[traj_j + 1]))
      }
      else{
        T_common <- Reduce(intersect, list(substr(traj_list[[traj_i]], 2, 2), substr(traj_list[[traj_j + 1]], 2, 2)))
        residual_cor[[length(residual_cor) + 1]] <- mxPath(from = paste0(traj_var[traj_i], T_common),
                                                           to = paste0(traj_var[traj_j + 1], T_common),
                                                           arrows = 2, free = T, values = starts[[3]][traj_i, traj_j + 1],
                                                           labels = paste0("residuals", traj_var[traj_i], traj_var[traj_j + 1]))
      }
      gf_cov_label[[traj_i + traj_j - 1]] <- matrix(paste0("psi", c("0s0s", "0s1s", "0s2s", "1s0s", "1s1s", "1s2s", "2s0s", "2s1s", "2s2s"),
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
  model_mx <- mxModel("Parallel Bilinear Spline Growth Model with Fixed Knots", type = "RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                             free = T, values = t(starts[[2]])[row(t(starts[[2]])) >= col(t(starts[[2]]))],
                             labels = t(multi_label)[row(t(multi_label)) >= col(t(multi_label))]),
                      outDef, outLoads1, outLoads2, outDef_L, outLoads1_L, outLoads2_L, loadings1, loadings2, loadings3,
                      gf_mean, gamma_mean, residual_var, residual_cor, func_L, grad_L, mean_s_L, psi_s_L, psi_btw_s_L, mean_L, psi_L, psi_btw_L)
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(original){
    model.para <- summary(model)$parameters[, c(1, 5, 6)]
    mean_est <- mean_se <- psi_est <- psi_se <- psi_btw_est <- psi_btw_se <- outcome_est <- outcome_se <- btw_est <- btw_se <- list()
    mug <- model.para$Estimate[c(grep("muknot", model.para$name))]
    mug_se <- model.para$Std.Error[c(grep("muknot", model.para$name))]
    for (traj in 1:length(traj_var)){
      mean_est[[length(mean_est) + 1]] <- c(mxEvalByName(paste0("mean", traj_var[traj]), model = model), mug[traj])
      mean_se[[length(mean_se) + 1]] <- c(mxSE(paste0("mean", traj_var[traj]), model, forceName = T), mug_se[traj])
      psi_est[[length(psi_est) + 1]] <- mxEvalByName(paste0("psi", traj_var[traj]), model = model)
      psi_se[[length(psi_se) + 1]] <- mxSE(paste0("psi", traj_var[traj]), model, forceName = T)
      outcome_est[[length(outcome_est) + 1]] <- c(unlist(mean_est[[traj]]), psi_est[[traj]][row(psi_est[[traj]]) >= col(psi_est[[traj]])])
      outcome_se[[length(outcome_se) + 1]] <- c(unlist(mean_se[[traj]]), psi_se[[traj]][row(psi_se[[traj]]) >= col(psi_se[[traj]])])
    }
    for (traj in 1:(length(traj_var) - 1)){
      psi_btw_est[[length(psi_btw_est) + 1]] <- mxEvalByName(paste0("psi", traj_var[traj], traj_var[traj + 1]), model = model)
      psi_btw_se[[length(psi_btw_se) + 1]] <- mxSE(paste0("psi", traj_var[traj], traj_var[traj + 1]), model, forceName = T)
      btw_est[[length(btw_est) + 1]] <- unlist(c(psi_btw_est))
      btw_se[[length(btw_se) + 1]] <- unlist(c(psi_btw_se))
    }
    model.est <- round(c(unlist(outcome_est), unlist(btw_est), model.para[grep("residuals", model.para$name), 2]), 4)
    model.se <- round(c(unlist(outcome_se), unlist(btw_se), model.para[grep("residuals", model.para$name), 3]), 4)
    estimate_out <- data.frame(Name = paraNames, Estimate = model.est, SE = model.se)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}


