getBLSGM_Med2 <- function(dat, M_records, Y_records, X_var, M_var, Y_var, t_var, res_ratio = rep(4, 2), btw_res = 0.3, diag_replace = F,
                          starts = NA, loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA){
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  if (any(is.na(starts))){
    #### Decide the initial values for predictor mean and variance
    mean_X <- mean(dat[, X_var], na.rm = T)
    var_X <- var(dat[, X_var], na.rm = T)
    starts.X <- list(mean_X, var_X)

    #### Decide the initial values for parameters of mediator trajectory
    M_BLSGM_F <- getBLSGM_Fixed(dat = dat, T_records = M_records, traj_var = M_var, t_var = t_var[1], res_ratio = res_ratio[1], original = F)
    M_gammaT <- M_BLSGM_F$mug$values
    M_time <- apply(dat[, paste0(t_var[1], M_records)], 2, mean)
    for (j in 1:length(M_time)){
      if (M_time[j] <= M_gammaT & M_time[j + 1] >= M_gammaT){
        M_gammaV <- apply(dat[, paste0(M_var, M_records)][, c(j, j + 1)], 1, mean)
        stop
      }
    }
    M_delta1 <- (M_gammaV - dat[, paste0(M_var, M_records[1])])/(rep(M_gammaT, nrow(dat)) - dat[, paste0(t_var[1], M_records[1])])
    M_delta2 <- (dat[, paste0(M_var, M_records[length(M_records)])] - M_gammaV)/(dat[, paste0(t_var[1], M_records[length(M_records)])]- rep(M_gammaT, nrow(dat)))
    M_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = M_delta1, X = dat[, X_var]), na.action = na.exclude)$coefficients)
    M_reg_r <- as.numeric(lm(M_gammaV ~ ., data = data.frame(M_gammaV = M_gammaV, X = dat[, X_var]),
                             na.action = na.exclude)$coefficients)
    M_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = M_delta2, X = dat[, X_var]), na.action = na.exclude)$coefficients)
    M_alpha0 <- c(M_reg_1[1], M_reg_r[1], M_reg_2[1], M_gammaT)
    beta_XtoM <- matrix(c(M_reg_1[2], M_reg_r[2], M_reg_2[2]), nrow = 3)
    M_growth_factor <- data.frame(M_eta1 = M_delta1, M_etar = M_gammaV, M_eta2 = M_delta2)
    M_psi0 <- cov(M_growth_factor)
    if (diag_replace){
      diag(M_psi0) <- c(M_BLSGM_F$psi$result[2, 2], M_BLSGM_F$psi_s$result[1, 1], M_BLSGM_F$psi$result[3, 3])
    }
    M_psi <- M_psi0  - beta_XtoM %*% var_X %*% t(beta_XtoM)
    starts.M <- list(M_alpha0, beta_XtoM, M_psi, M_BLSGM_F$S$values[1, 1])

    #### Decide the initial values for parameters of outcome trajectory
    Y_BLSGM_F <- getBLSGM_Fixed(dat = dat, T_records = Y_records, traj_var = Y_var, t_var = t_var[2], res_ratio = res_ratio[2], original = F)
    Y_gammaT <- Y_BLSGM_F$mug$values
    Y_time <- apply(dat[, paste0(t_var[2], Y_records)], 2, mean)
    for (j in 1:length(Y_time)){
      if (Y_time[j] <= Y_gammaT & Y_time[j + 1] >= Y_gammaT){
        Y_gammaV <- apply(dat[, paste0(Y_var, Y_records)][, c(j, j + 1)], 1, mean)
        stop
      }
    }
    Y_delta1 <- (Y_gammaV - dat[, paste0(Y_var, Y_records[1])])/(rep(Y_gammaT, nrow(dat)) - dat[, paste0(t_var[2], Y_records[1])])
    Y_delta2 <- (dat[, paste0(Y_var, Y_records[length(Y_records)])] - Y_gammaV)/(dat[, paste0(t_var[2], Y_records[length(Y_records)])]- rep(Y_gammaT, nrow(dat)))
    Y_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_delta1, X = dat[, X_var], M_slp1 = M_delta1),
                             na.action = na.exclude)$coefficients)
    Y_reg_r <- as.numeric(lm(Y_gammaV ~ ., data = data.frame(Y_gammaV = Y_gammaV, X = dat[, X_var], M_slp1 = M_delta1, M_gammaV = M_gammaV),
                             na.action = na.exclude)$coefficients)
    Y_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_delta2, X = dat[, X_var],
                                                               M_slp1 = M_delta1, M_gammaV = M_gammaV, M_slp2 = M_delta2),
                             na.action = na.exclude)$coefficients)
    Y_alpha0 <- c(Y_reg_1[1], Y_reg_r[1], Y_reg_2[1], Y_gammaT)
    beta_XtoY <- matrix(c(Y_reg_1[2], Y_reg_r[2], Y_reg_2[2]), nrow = 3)
    beta_MtoY <- matrix(c(Y_reg_1[3], 0, 0, Y_reg_r[3:4], 0, Y_reg_2[3:5]), byrow = T, nrow = 3, ncol = 3)
    Y_growth_factor <- data.frame(Y_eta1 = Y_delta1, Y_etar = Y_gammaV, Y_eta2 = Y_delta2)
    Y_psi0 <- cov(Y_growth_factor)
    if (diag_replace){
      diag(Y_psi0) <- c(Y_BLSGM_F$psi$result[2, 2], Y_BLSGM_F$psi_s$result[1, 1], Y_BLSGM_F$psi$result[3, 3])
    }
    Y_psi <- Y_psi0 - beta_MtoY %*% M_psi %*% t(beta_MtoY) - beta_XtoY %*% var_X %*% t(beta_XtoY)
    starts.Y <- list(Y_alpha0, beta_XtoY, beta_MtoY, Y_psi, Y_BLSGM_F$S$values[1, 1])
    starts <- list(starts.X, starts.M, starts.Y)
  }
  ### Define manifest variables
  traj_var <- c(Y_var, M_var)
  T_records <- list(Y_records, M_records)
  traj_list <- list()
  for (traj in 1:length(traj_var)){
    traj_list[[length(traj_list) + 1]] <- paste0(traj_var[traj], T_records[[traj]])
  }
  manifests <- c(unlist(traj_list), X_var)
  ### Define latent variables
  latents <- c("etaY1", "etaYr", "etaY2", "etaM1", "etaMr", "etaM2")
  outDefM <- outDefY <- list(); outLoadsY1 <- outLoadsY2 <- outLoadsM1 <- outLoadsM2 <- list()
  for (j in M_records){
    outDefM[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[1], j), name = paste0("t", M_var, j))
    outLoadsM1[[j]] <- mxAlgebraFromString(paste0("min(0, t", M_var, j, " - mugM)"), name = paste0("L1", j, "M"))
    outLoadsM2[[j]] <- mxAlgebraFromString(paste0("max(0, t", M_var, j, " - mugM)"), name = paste0("L2", j, "M"))
  }
  for (j in Y_records){
    outDefY[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[2], j), name = paste0("t", Y_var, j))
    outLoadsY1[[j]] <- mxAlgebraFromString(paste0("min(0, t", Y_var, j, " - mugY)"), name = paste0("L1", j, "Y"))
    outLoadsY2[[j]] <- mxAlgebraFromString(paste0("max(0, t", Y_var, j, " - mugY)"), name = paste0("L2", j, "Y"))
  }
  residual_cor <- list()
  for (traj_i in 1:(length(traj_var) - 1)){
    for (traj_j in traj_i:(length(traj_var) - 1)){
      #### Define the covariances of residuals
      if (setequal(readr::parse_number(traj_list[[traj_i]]), readr::parse_number(traj_list[[traj_j + 1]]))){
        residual_cor[[length(residual_cor) + 1]] <- mxPath(from = traj_list[[traj_i]], to = traj_list[[traj_j + 1]],
                                                           arrows = 2, free = T, values = btw_res[traj_i + traj_j - 1] * sqrt(starts[[3]][[5]] * starts[[2]][[4]]),
                                                           labels = paste0("residuals", traj_var[traj_i], traj_var[traj_j + 1]))
      }
      else{
        T_common <- Reduce(intersect, list(readr::parse_number(traj_list[[traj_i]]), readr::parse_number(traj_list[[traj_j + 1]])))
        residual_cor[[length(residual_cor) + 1]] <- mxPath(from = paste0(traj_var[traj_i], T_common),
                                                           to = paste0(traj_var[traj_j + 1], T_common),
                                                           arrows = 2, free = T, values = btw_res[traj_i + traj_j - 1] * sqrt(starts[[3]][[5]] * starts[[2]][[4]]),
                                                           labels = paste0("residuals", traj_var[traj_i], traj_var[traj_j + 1]))
      }
    }
  }
  ### Create a mxModel object
  model_mx <- mxModel("Mediation Process in Bivariate Nonlinear Growth Model with Fixed Knots", type = "RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define factor loadings from latent variables to manifests
                      mxPath(from = "etaY1", to = paste0(Y_var, Y_records), arrows = 1, free = F, values = 0,
                             labels = paste0("L1", Y_records, "Y[1,1]")),
                      mxPath(from = "etaYr", to = paste0(Y_var, Y_records), arrows = 1, free = F, values = 1),
                      mxPath(from = "etaY2", to = paste0(Y_var, Y_records), arrows = 1, free = F, values = 0,
                             labels = paste0("L2", Y_records, "Y[1,1]")),
                      mxPath(from = "etaM1", to = paste0(M_var, M_records), arrows = 1, free = F, values = 0,
                             labels = paste0("L1", M_records, "M[1,1]")),
                      mxPath(from = "etaMr", to = paste0(M_var, M_records), arrows = 1, free = F, values = 1),
                      mxPath(from = "etaM2", to = paste0(M_var, M_records), arrows = 1, free = F, values = 0,
                             labels = paste0("L2", M_records, "M[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = paste0(Y_var, Y_records), to = paste0(Y_var, Y_records),
                             arrows = 2, free = T, values = starts[[3]][[5]],
                             labels = "residualsY"),
                      mxPath(from = paste0(M_var, M_records), to = paste0(M_var, M_records),
                             arrows = 2, free = T, values = starts[[2]][[4]],
                             labels = "residualsM"),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents, arrows = 1, free = T, values = c(starts[[3]][[1]][1:3], starts[[2]][[1]][1:3]),
                             labels = c("alphaY1", "alphaYr", "alphaY2", "alphaM1", "alphaMr", "alphaM2")),
                      #### Define var-cov matrix of within-construct growth factors
                      mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[[3]][[4]][row(starts[[3]][[4]]) >= col(starts[[3]][[4]])],
                             labels = c("psiY1Y1", "psiY1Yr", "psiY1Y2", "psiYrYr", "psiYrY2", "psiY2Y2")),
                      mxPath(from = latents[4:6], to = latents[4:6], arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[[2]][[3]][row(starts[[2]][[3]]) >= col(starts[[2]][[3]])],
                             labels = c("psiM1M1", "psiM1Mr", "psiM1M2", "psiMrMr", "psiMrM2", "psiM2M2")),
                      #### Include time-invariant covariates
                      ##### Means
                      mxPath(from = "one", to = "X", arrows = 1, free = T, values = starts[[1]][[1]], labels = c("muX")),
                      mxPath(from = "X", to = "X", connect = "unique.pairs", arrows = 2, free = T, values = c(starts[[1]][[2]]), labels = "phi11"),
                      ##### Regression coefficients (Direct effect on growth factors of X)
                      mxPath(from = "X", to = latents, arrows = 1, free = T, values = c(starts[[3]][[2]], starts[[2]][[2]]),
                             labels = paste0("beta", rep(c("Y", "M"), each = 3), rep(c(1, "r", 2), 2))),
                      #### Define coefficients that contribute to indirect effects between-construct growth factors
                      ##### Slope 1 of M to slope 1 of Y
                      mxPath(from = latents[4], to = latents[1], arrows = 1, free = T, values = starts[[3]][[3]][1, 1],
                             labels = "betaM1Y1"),
                      ##### Slope 1 of M to intercept of Y
                      mxPath(from = latents[4], to = latents[2], arrows = 1, free = T, values = starts[[3]][[3]][2, 1],
                             labels = "betaM1Yr"),
                      ##### Slope 1 of M to slope 2 of Y
                      mxPath(from = latents[4], to = latents[3], arrows = 1, free = T, values = starts[[3]][[3]][3, 1],
                             labels = "betaM1Y2"),
                      ##### Intercept of M to intercept of Y
                      mxPath(from = latents[5], to = latents[2], arrows = 1, free = T, values = starts[[3]][[3]][2, 2],
                             labels = "betaMrYr"),
                      ##### Intercept of M to slope 2 of Y
                      mxPath(from = latents[5], to = latents[3], arrows = 1, free = T, values = starts[[3]][[3]][3, 2],
                             labels = "betaMrY2"),
                      ##### Slope 2 of M to slope 2 of Y
                      mxPath(from = latents[6], to = latents[3], arrows = 1, free = T, values = starts[[3]][[3]][3, 3],
                             labels = "betaM2Y2"),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = starts[[3]][[1]][4], labels = "muknot_Y", name = "mugY"),
                      mxMatrix("Full", 1, 1, free = T, values = starts[[2]][[1]][4], labels = "muknot_M", name = "mugM"),
                      outDefM, outDefY, outLoadsY1, outLoadsY2, outLoadsM1, outLoadsM2, residual_cor,

                      #### Calculate the mean vector of outcome Y
                      ##### Intercept coefficients of mediator M
                      mxAlgebra(rbind(alphaM1, alphaMr, alphaM2), name = "alphaM"),
                      ##### Intercept coefficients of outcome Y
                      mxAlgebra(rbind(alphaY1, alphaYr, alphaY2), name = "alphaY"),
                      ##### Coefficients form X to mediator M
                      mxAlgebra(rbind(betaM1, betaMr, betaM2), name = "beta_xm"),
                      ##### Coefficients form X to mediator Y
                      mxAlgebra(rbind(betaY1, betaYr, betaY2), name = "beta_xy"),
                      ##### Coefficients from mediator M to outcome Y
                      mxAlgebra(rbind(cbind(betaM1Y1, 0, 0),
                                      cbind(betaM1Yr, betaMrYr, 0),
                                      cbind(betaM1Y2, betaMrY2, betaM2Y2)), name = "beta_my"),
                      mxAlgebra(alphaM + beta_xm %*% muX, name = "muetaM"),
                      mxAlgebra(alphaY + beta_my %*% muetaM + beta_xy %*% muX, name = "muetaY"),
                      mxAlgebra(betaM1Y1 * betaM1, name = "mediator_11"),
                      mxAlgebra(betaM1Yr * betaM1, name = "mediator_1r"),
                      mxAlgebra(betaM1Y2 * betaM1, name = "mediator_12"),
                      mxAlgebra(betaMrYr * betaMr, name = "mediator_rr"),
                      mxAlgebra(betaMrY2 * betaMr, name = "mediator_r2"),
                      mxAlgebra(betaM2Y2 * betaM2, name = "mediator_22"),
                      mxAlgebra(beta_my %*% beta_xm + beta_xy, name = "total"))
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(original){
    model.para <- summary(model)$parameters[, c(1, 5, 6)]
    model.est <- round(c(model$muetaY$result, model.para[model.para$name == "muknot_Y", 2],
                         model$muetaM$result, model.para[model.para$name == "muknot_M", 2],
                         model.para[grep("psiY1", model.para$name), 2],
                         model.para[grep("psiYr", model.para$name), 2],
                         model.para[grep("psiY2", model.para$name), 2],
                         model.para[grep("psiM1", model.para$name), 2],
                         model.para[grep("psiMr", model.para$name), 2],
                         model.para[grep("psiM2", model.para$name), 2],
                         model.para[model.para$name == "muX", 2],
                         model.para[model.para$name == "phi11", 2],
                         model$mediator_11$result, model$mediator_1r$result, model$mediator_rr$result,
                         model$mediator_12$result, model$mediator_r2$result, model$mediator_22$result,
                         model$total$result,
                         model.para[c(1:12, 13, 15, 14), 2]), 4)

    model.se <- round(c(mxSE(muetaY, model), model.para[model.para$name == "muknot_Y", 3],
                        mxSE(muetaM, model), model.para[model.para$name == "muknot_M", 3],
                        model.para[grep("psiY1", model.para$name), 3],
                        model.para[grep("psiYr", model.para$name), 3],
                        model.para[grep("psiY2", model.para$name), 3],
                        model.para[grep("psiM1", model.para$name), 3],
                        model.para[grep("psiMr", model.para$name), 3],
                        model.para[grep("psiM2", model.para$name), 3],
                        model.para[model.para$name == "muX", 3],
                        model.para[model.para$name == "phi11", 3],
                        mxSE(mediator_11, model), mxSE(mediator_1r, model), mxSE(mediator_rr, model),
                        mxSE(mediator_12, model), mxSE(mediator_r2, model), mxSE(mediator_22, model),
                        mxSE(total, model),
                        model.para[c(1:12, 13, 15, 14), 3]), 4)
    estimate_out <- data.frame(Name = paraNames, Estimate = model.est, Std.Error = model.se)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}


