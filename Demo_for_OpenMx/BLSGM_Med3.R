getBLSGM_Med3 <- function(dat, X_records, M_records, Y_records, X_var, M_var, Y_var, t_var, res_ratio = rep(4, 3), btw_res = rep(0.3, 3),
                          diag_replace = F, starts = NA, loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA){
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  if (any(is.na(starts))){
    #### Decide the initial values for predictor mean and variance
    X_BLSGM_F <- getBLSGM_Fixed(dat = dat, T_records = X_records, traj_var = X_var, t_var = t_var[1], res_ratio = res_ratio[1], original = F)
    X_gammaT <- X_BLSGM_F$mug$values
    X_mean0 <- c(X_BLSGM_F$mean$result[2], X_BLSGM_F$mean_s$result[1], X_BLSGM_F$mean$result[3], X_gammaT)
    X_time <- apply(dat[, paste0(t_var[1], X_records)], 2, mean)
    for (j in 1:length(X_time)){
      if (X_time[j] <= X_gammaT & X_time[j + 1] >= X_gammaT){
        X_gammaV <- apply(dat[, paste0(X_var, X_records)][, c(j, j + 1)], 1, mean)
        stop
      }
    }
    X_delta1 <- (X_gammaV - dat[, paste0(X_var, X_records[1])])/(rep(X_gammaT, nrow(dat)) - dat[, paste0(t_var[1], X_records[1])])
    X_delta2 <- (dat[, paste0(X_var, X_records[length(X_records)])] - X_gammaV)/(dat[, paste0(t_var[1], X_records[length(X_records)])]- rep(X_gammaT, nrow(dat)))
    X_growth_factor <- data.frame(X_eta1 = X_delta1, X_etar = X_gammaV, X_eta2 = X_delta2)
    X_psi0 <- cov(X_growth_factor)
    if (diag_replace){
      diag(X_psi0) <- c(X_BLSGM_F$psi$result[2, 2], X_BLSGM_F$psi_s$result[1, 1], X_BLSGM_F$psi$result[3, 3])
    }
    starts.X <- list(X_mean0, X_psi0, X_BLSGM_F$S$values[1, 1])

    #### Decide the initial values for parameters of mediator trajectory
    M_BLSGM_F <- getBLSGM_Fixed(dat = dat, T_records = M_records, traj_var = M_var, t_var = t_var[2], res_ratio = res_ratio[2], original = F)
    M_gammaT <- M_BLSGM_F$mug$values
    M_time <- apply(dat[, paste0(t_var[2], M_records)], 2, mean)
    for (j in 1:length(M_time)){
      if (M_time[j] <= M_gammaT & M_time[j + 1] >= M_gammaT){
        M_gammaV <- apply(dat[, paste0(M_var, M_records)][, c(j, j + 1)], 1, mean)
        stop
      }
    }
    M_delta1 <- (M_gammaV - dat[, paste0(M_var, M_records[1])])/(rep(M_gammaT, nrow(dat)) - dat[, paste0(t_var[2], M_records[1])])
    M_delta2 <- (dat[, paste0(M_var, M_records[length(M_records)])] - M_gammaV)/(dat[, paste0(t_var[2], M_records[length(M_records)])]- rep(M_gammaT, nrow(dat)))
    M_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = M_delta1, X_slp1 = X_delta1), na.action = na.exclude)$coefficients)
    M_reg_r <- as.numeric(lm(M_gammaV ~ ., data = data.frame(M_gammaV = M_gammaV, X_slp1 = X_delta1, X_gammaV = X_gammaV),
                             na.action = na.exclude)$coefficients)
    M_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = M_delta2, X_slp1 = X_delta1, X_gammaV = X_gammaV,
                                                               X_slp2 = X_delta2), na.action = na.exclude)$coefficients)
    M_alpha0 <- c(M_reg_1[1], M_reg_r[1], M_reg_2[1], M_gammaT)
    beta_XtoM <- matrix(c(M_reg_1[-1], 0, 0, M_reg_r[-1], 0, M_reg_2[-1]), byrow = T, nrow = 3, ncol = 3)
    M_growth_factor <- data.frame(M_eta1 = M_delta1, M_etar = M_gammaV, M_eta2 = M_delta2)
    M_psi0 <- cov(M_growth_factor)
    if (diag_replace){
      diag(M_psi0) <- c(M_BLSGM_F$psi$result[2, 2], M_BLSGM_F$psi_s$result[1, 1], M_BLSGM_F$psi$result[3, 3])
    }
    M_psi <- M_psi0  - beta_XtoM %*% X_psi0 %*% t(beta_XtoM)
    starts.M <- list(M_alpha0, beta_XtoM, M_psi, M_BLSGM_F$S$values[1, 1])

    #### Decide the initial values for parameters of outcome trajectory
    Y_BLSGM_F <- getBLSGM_Fixed(dat = dat, T_records = Y_records, traj_var = Y_var, t_var = t_var[3], res_ratio = res_ratio[3], original = F)
    Y_gammaT <- Y_BLSGM_F$mug$values
    Y_time <- apply(dat[, paste0(t_var[3], Y_records)], 2, mean)
    for (j in 1:length(Y_time)){
      if (Y_time[j] <= Y_gammaT & Y_time[j + 1] >= Y_gammaT){
        Y_gammaV <- apply(dat[, paste0(Y_var, Y_records)][, c(j, j + 1)], 1, mean)
        stop
      }
    }
    Y_delta1 <- (Y_gammaV - dat[, paste0(Y_var, Y_records[1])])/(rep(Y_gammaT, nrow(dat)) - dat[, paste0(t_var[3], Y_records[1])])
    Y_delta2 <- (dat[, paste0(Y_var, Y_records[length(Y_records)])] - Y_gammaV)/(dat[, paste0(t_var[3], Y_records[length(Y_records)])]- rep(Y_gammaT, nrow(dat)))
    Y_reg_1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_delta1, X_slp1 = X_delta1, M_slp1 = M_delta1),
                             na.action = na.exclude)$coefficients)
    Y_reg_r <- as.numeric(lm(Y_gammaV ~ ., data = data.frame(Y_gammaV = Y_gammaV, X_slp1 = X_delta1, X_gammaV = X_gammaV,
                                                             M_slp1 = M_delta1, M_gammaV = M_gammaV),
                             na.action = na.exclude)$coefficients)
    Y_reg_2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = Y_delta2, X_slp1 = X_delta1, X_gammaV = X_gammaV, X_slp2 = X_delta2,
                                                               M_slp1 = M_delta1, M_gammaV = M_gammaV, M_slp2 = M_delta2),
                             na.action = na.exclude)$coefficients)
    Y_alpha0 <- c(Y_reg_1[1], Y_reg_r[1], Y_reg_2[1], Y_gammaT)
    beta_XtoY <- matrix(c(Y_reg_1[2], 0, 0, Y_reg_r[2:3], 0, Y_reg_2[2:4]), byrow = T, nrow = 3, ncol = 3)
    beta_MtoY <- matrix(c(Y_reg_1[3], 0, 0, Y_reg_r[4:5], 0, Y_reg_2[5:7]), byrow = T, nrow = 3, ncol = 3)
    Y_growth_factor <- data.frame(Y_eta1 = Y_delta1, Y_etar = Y_gammaV, Y_eta2 = Y_delta2)
    Y_psi0 <- cov(Y_growth_factor)
    if (diag_replace){
      diag(Y_psi0) <- c(Y_BLSGM_F$psi$result[2, 2], Y_BLSGM_F$psi_s$result[1, 1], Y_BLSGM_F$psi$result[3, 3])
    }
    Y_psi <- Y_psi0 - beta_MtoY %*% M_psi %*% t(beta_MtoY) - beta_XtoY %*% X_psi0 %*% t(beta_XtoY)
    starts.Y <- list(Y_alpha0, beta_XtoY, beta_MtoY, Y_psi, Y_BLSGM_F$S$values[1, 1])
    starts <- list(starts.X, starts.M, starts.Y)
  }
  ### Define manifest variables
  traj_var <- c(Y_var, M_var, X_var)
  T_records <- list(Y_records, M_records, X_records)
  traj_list <- list()
  for (traj in 1:length(traj_var)){
    traj_list[[length(traj_list) + 1]] <- paste0(traj_var[traj], T_records[[traj]])
  }
  manifests <- unlist(traj_list)
  ### Define latent variables: growth factors for Reading and Math IRT trajectories
  latents <- c("etaY1", "etaYr", "etaY2", "etaM1", "etaMr", "etaM2", "etaX1", "etaXr", "etaX2")
  outDefX <- outDefM <- outDefY <- list(); outLoadsY1 <- outLoadsY2 <- outLoadsM1 <- outLoadsM2 <- outLoadsX1 <- outLoadsX2 <- list()
  for (j in X_records){
    outDefX[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[1], j), name = paste0("t", X_var, j))
    outLoadsX1[[j]] <- mxAlgebraFromString(paste0("min(0, t", X_var, j, " - mugX)"), name = paste0("L1", j, "X"))
    outLoadsX2[[j]] <- mxAlgebraFromString(paste0("max(0, t", X_var, j, " - mugX)"), name = paste0("L2", j, "X"))
  }
  for (j in M_records){
    outDefM[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[2], j), name = paste0("t", M_var, j))
    outLoadsM1[[j]] <- mxAlgebraFromString(paste0("min(0, t", M_var, j, " - mugM)"), name = paste0("L1", j, "M"))
    outLoadsM2[[j]] <- mxAlgebraFromString(paste0("max(0, t", M_var, j, " - mugM)"), name = paste0("L2", j, "M"))
  }
  for (j in Y_records){
    outDefY[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[3], j), name = paste0("t", Y_var, j))
    outLoadsY1[[j]] <- mxAlgebraFromString(paste0("min(0, t", Y_var, j, " - mugY)"), name = paste0("L1", j, "Y"))
    outLoadsY2[[j]] <- mxAlgebraFromString(paste0("max(0, t", Y_var, j, " - mugY)"), name = paste0("L2", j, "Y"))
  }
  var_1 <- c(starts[[3]][[5]], starts[[3]][[5]], starts[[2]][[4]])
  var_2 <- c(starts[[2]][[4]], starts[[1]][[3]], starts[[1]][[3]])
  residual_cor <- list()
  for (traj_i in 1:(length(traj_var) - 1)){
    for (traj_j in traj_i:(length(traj_var) - 1)){
      #### Define the covariances of residuals
      if (setequal(readr::parse_number(traj_list[[traj_i]]), readr::parse_number(traj_list[[traj_j + 1]]))){
        residual_cor[[length(residual_cor) + 1]] <- mxPath(from = traj_list[[traj_i]], to = traj_list[[traj_j + 1]],
                                                           arrows = 2, free = T, values = btw_res[traj_i + traj_j - 1] * sqrt(var_1[traj_i] * var_2[traj_j]),
                                                           labels = paste0("residuals", traj_var[traj_i], traj_var[traj_j + 1]))
      }
      else{
        T_common <- Reduce(intersect, list(readr::parse_number(traj_list[[traj_i]]), readr::parse_number(traj_list[[traj_j + 1]])))
        residual_cor[[length(residual_cor) + 1]] <- mxPath(from = paste0(traj_var[traj_i], T_common),
                                                           to = paste0(traj_var[traj_j + 1], T_common),
                                                           arrows = 2, free = T, values = btw_res[traj_i + traj_j - 1] * sqrt(var_1[traj_i] * var_2[traj_j]),
                                                           labels = paste0("residuals", traj_var[traj_i], traj_var[traj_j + 1]))
      }
    }
  }

  ### Create a mxModel object
  model_mx <- mxModel("Mediation Process in Trivariate Nonlinear Growth Model with Fixed Knots", type = "RAM",
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
                      mxPath(from = "etaX1", to = paste0(X_var, X_records), arrows = 1, free = F, values = 0,
                             labels = paste0("L1", X_records, "X[1,1]")),
                      mxPath(from = "etaXr", to = paste0(X_var, X_records), arrows = 1, free = F, values = 1),
                      mxPath(from = "etaX2", to = paste0(X_var, X_records), arrows = 1, free = F, values = 0,
                             labels = paste0("L2", X_records, "X[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = paste0(Y_var, Y_records), to = paste0(Y_var, Y_records),
                             arrows = 2, free = T, values = starts[[3]][[5]], labels = "residualsY"),
                      mxPath(from = paste0(M_var, M_records), to = paste0(M_var, M_records),
                             arrows = 2, free = T, values = starts[[2]][[4]], labels = "residualsM"),
                      mxPath(from = paste0(X_var, X_records), to = paste0(X_var, X_records),
                             arrows = 2, free = T, values = starts[[1]][[3]], labels = "residualsX"),

                      #### Define means of latent variables
                      mxPath(from = "one", to = latents, arrows = 1, free = T,
                             values = c(starts[[3]][[1]][1:3], starts[[2]][[1]][1:3], starts[[1]][[1]][1:3]),
                             labels = c("alphaY1", "alphaYr", "alphaY2",
                                        "alphaM1", "alphaMr", "alphaM2",
                                        "alphaX1", "alphaXr", "alphaX2")),
                      #### Define var-cov matrix of within-construct growth factors
                      mxPath(from = latents[1:3], to = latents[1:3], arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[[3]][[4]][row(starts[[3]][[4]]) >= col(starts[[3]][[4]])],
                             labels = c("psiY1Y1", "psiY1Yr", "psiY1Y2",
                                        "psiYrYr", "psiYrY2",
                                        "psiY2Y2")),
                      mxPath(from = latents[4:6], to = latents[4:6], arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[[2]][[3]][row(starts[[2]][[3]]) >= col(starts[[2]][[3]])],
                             labels = c("psiM1M1", "psiM1Mr", "psiM1M2",
                                        "psiMrMr", "psiMrM2",
                                        "psiM2M2")),
                      mxPath(from = latents[7:9], to = latents[7:9], arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[[1]][[2]][row(starts[[1]][[2]]) >= col(starts[[1]][[2]])],
                             labels = c("psiX1X1", "psiX1Xr", "psiX1X2",
                                        "psiXrXr", "psiXrX2",
                                        "psiX2X2")),

                      #### Define coefficients that contribute to indirect effects between-construct growth factors
                      ##### Slope 1 of X to slope 1 of Y
                      mxPath(from = latents[7], to = latents[1], arrows = 1, free = T, values = starts[[3]][[2]][1, 1],
                             labels = "betaX1Y1"),
                      ##### Slope 1 of X to intercept of Y
                      mxPath(from = latents[7], to = latents[2], arrows = 1, free = T, values = starts[[3]][[2]][2, 1],
                             labels = "betaX1Yr"),
                      ##### Slope 1 of X to slope 2 of Y
                      mxPath(from = latents[7], to = latents[3], arrows = 1, free = T, values = starts[[3]][[2]][3, 1],
                             labels = "betaX1Y2"),
                      ##### Intercept of X to intercept of Y
                      mxPath(from = latents[8], to = latents[2], arrows = 1, free = T, values = starts[[3]][[2]][2, 2],
                             labels = "betaXrYr"),
                      ##### Intercept of X to slope 2 of Y
                      mxPath(from = latents[8], to = latents[3], arrows = 1, free = T, values = starts[[3]][[2]][3, 2],
                             labels = "betaXrY2"),
                      ##### Slope 2 of X to slope 2 of Y
                      mxPath(from = latents[9], to = latents[3], arrows = 1, free = T, values = starts[[3]][[2]][3, 3],
                             labels = "betaX2Y2"),

                      ##### Slope 1 of X to slope 1 of M
                      mxPath(from = latents[7], to = latents[4], arrows = 1, free = T, values = starts[[2]][[2]][1, 1],
                             labels = "betaX1M1"),
                      ##### Slope 1 of X to intercept of M
                      mxPath(from = latents[7], to = latents[5], arrows = 1, free = T, values = starts[[2]][[2]][2, 1],
                             labels = "betaX1Mr"),
                      ##### Slope 1 of X to slope 2 of M
                      mxPath(from = latents[7], to = latents[6], arrows = 1, free = T, values = starts[[2]][[2]][3, 1],
                             labels = "betaX1M2"),
                      ##### Intercept of X to intercept of M
                      mxPath(from = latents[8], to = latents[5], arrows = 1, free = T, values = starts[[2]][[2]][2, 2],
                             labels = "betaXrMr"),
                      ##### Intercept of X to slope 2 of M
                      mxPath(from = latents[8], to = latents[6], arrows = 1, free = T, values = starts[[2]][[2]][3, 2],
                             labels = "betaXrM2"),
                      ##### Slope 2 of X to slope 2 of M
                      mxPath(from = latents[9], to = latents[6], arrows = 1, free = T, values = starts[[2]][[2]][3, 3],
                             labels = "betaX2M2"),

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
                      mxMatrix("Full", 1, 1, free = T, values = starts[[3]][[1]][4],
                               labels = "muknot_Y", name = "mugY"),
                      mxMatrix("Full", 1, 1, free = T, values = starts[[2]][[1]][4],
                               labels = "muknot_M", name = "mugM"),
                      mxMatrix("Full", 1, 1, free = T, values = starts[[1]][[1]][4],
                               labels = "muknot_X", name = "mugX"),
                      outDefY, outDefM, outDefX, outLoadsY1, outLoadsY2, outLoadsM1, outLoadsM2, outLoadsX1, outLoadsX2, residual_cor,

                      #### Calculate the mean vector of outcome Y
                      ##### Intercept coefficients of intervention X
                      mxAlgebra(rbind(alphaX1, alphaXr, alphaX2), name = "alphaX"),
                      ##### Intercept coefficients of mediator M
                      mxAlgebra(rbind(alphaM1, alphaMr, alphaM2), name = "alphaM"),
                      ##### Intercept coefficients of outcome Y
                      mxAlgebra(rbind(alphaY1, alphaYr, alphaY2), name = "alphaY"),
                      ##### Coefficients form X to mediator M
                      mxAlgebra(rbind(cbind(betaX1M1, 0, 0),
                                      cbind(betaX1Mr, betaXrMr, 0),
                                      cbind(betaX1M2, betaXrM2, betaX2M2)), name = "beta_xm"),
                      ##### Coefficients form X to mediator Y
                      mxAlgebra(rbind(cbind(betaX1Y1, 0, 0),
                                      cbind(betaX1Yr, betaXrYr, 0),
                                      cbind(betaX1Y2, betaXrY2, betaX2Y2)), name = "beta_xy"),
                      ##### Coefficients from mediator M to outcome Y
                      mxAlgebra(rbind(cbind(betaM1Y1, 0, 0),
                                      cbind(betaM1Yr, betaMrYr, 0),
                                      cbind(betaM1Y2, betaMrY2, betaM2Y2)), name = "beta_my"),
                      mxAlgebra(alphaM + beta_xm %*% alphaX, name = "muetaM"),
                      mxAlgebra(alphaY + beta_my %*% muetaM + beta_xy %*% alphaX, name = "muetaY"),
                      ##### Inference of indirect effects
                      ###### X1--M1--Y1
                      mxAlgebra(betaX1M1 * betaM1Y1, name = "mediator_111"),
                      ###### Xr--Mr--Yr
                      mxAlgebra(betaXrMr * betaMrYr, name = "mediator_rrr"),
                      ###### X2--M2--Y2
                      mxAlgebra(betaX2M2 * betaM2Y2, name = "mediator_222"),
                      ###### X1--M1--Yr
                      mxAlgebra(betaX1M1 * betaM1Yr, name = "mediator_11r"),
                      ###### X1--Mr--Yr
                      mxAlgebra(betaX1Mr * betaMrYr, name = "mediator_1rr"),
                      ###### X1--M1--Y2
                      mxAlgebra(betaX1M1 * betaM1Y2, name = "mediator_112"),
                      ###### X1--Mr--Y2
                      mxAlgebra(betaX1Mr * betaMrY2, name = "mediator_1r2"),
                      ###### X1--M2--Y2
                      mxAlgebra(betaX1M2 * betaM2Y2, name = "mediator_122"),
                      ###### Xr--Mr--Y2
                      mxAlgebra(betaXrMr * betaMrY2, name = "mediator_rr2"),
                      ###### Xr--M2--Y2
                      mxAlgebra(betaXrM2 * betaM2Y2, name = "mediator_r22"),
                      ##### Total effect
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
                         model$alphaX$result, model.para[model.para$name == "muknot_X", 2],
                         model.para[model.para$name == "psiY1Y1", 2],
                         model.para[model.para$name == "psiY1Yr", 2],
                         model.para[model.para$name == "psiY1Y2", 2],
                         model.para[model.para$name == "psiYrYr", 2],
                         model.para[model.para$name == "psiYrY2", 2],
                         model.para[model.para$name == "psiY2Y2", 2],
                         model.para[model.para$name == "psiM1M1", 2],
                         model.para[model.para$name == "psiM1Mr", 2],
                         model.para[model.para$name == "psiM1M2", 2],
                         model.para[model.para$name == "psiMrMr", 2],
                         model.para[model.para$name == "psiMrM2", 2],
                         model.para[model.para$name == "psiM2M2", 2],
                         model.para[model.para$name == "psiX1X1", 2],
                         model.para[model.para$name == "psiX1Xr", 2],
                         model.para[model.para$name == "psiX1X2", 2],
                         model.para[model.para$name == "psiXrXr", 2],
                         model.para[model.para$name == "psiXrX2", 2],
                         model.para[model.para$name == "psiX2X2", 2],

                         model.para[model.para$name == "betaX1Y1", 2],
                         model.para[model.para$name == "betaX1Yr", 2],
                         model.para[model.para$name == "betaX1Y2", 2],
                         model.para[model.para$name == "betaXrYr", 2],
                         model.para[model.para$name == "betaXrY2", 2],
                         model.para[model.para$name == "betaX2Y2", 2],
                         model.para[model.para$name == "betaX1M1", 2],
                         model.para[model.para$name == "betaX1Mr", 2],
                         model.para[model.para$name == "betaX1M2", 2],
                         model.para[model.para$name == "betaXrMr", 2],
                         model.para[model.para$name == "betaXrM2", 2],
                         model.para[model.para$name == "betaX2M2", 2],
                         model.para[model.para$name == "betaM1Y1", 2],
                         model.para[model.para$name == "betaM1Yr", 2],
                         model.para[model.para$name == "betaM1Y2", 2],
                         model.para[model.para$name == "betaMrYr", 2],
                         model.para[model.para$name == "betaMrY2", 2],
                         model.para[model.para$name == "betaM2Y2", 2],
                         model$mediator_111$result, model$mediator_rrr$result, model$mediator_222$result,
                         model$mediator_11r$result, model$mediator_1rr$result,
                         model$mediator_112$result, model$mediator_1r2$result, model$mediator_122$result,
                         model$mediator_rr2$result, model$mediator_r22$result,
                         c(model$total$result[1:3, 1], model$total$result[2:3, 2], model$total$result[3, 3]),
                         model.para[c(19, 21, 24, 20, 22, 23), 2]), 4)

    model.se <- round(c(mxSE(muetaY, model), model.para[model.para$name == "muknot_Y", 3],
                        mxSE(muetaM, model), model.para[model.para$name == "muknot_M", 3],
                        mxSE(alphaX, model), model.para[model.para$name == "muknot_X", 3],
                        model.para[model.para$name == "psiY1Y1", 3],
                        model.para[model.para$name == "psiY1Yr", 3],
                        model.para[model.para$name == "psiY1Y2", 3],
                        model.para[model.para$name == "psiYrYr", 3],
                        model.para[model.para$name == "psiYrY2", 3],
                        model.para[model.para$name == "psiY2Y2", 3],
                        model.para[model.para$name == "psiM1M1", 3],
                        model.para[model.para$name == "psiM1Mr", 3],
                        model.para[model.para$name == "psiM1M2", 3],
                        model.para[model.para$name == "psiMrMr", 3],
                        model.para[model.para$name == "psiMrM2", 3],
                        model.para[model.para$name == "psiM2M2", 3],
                        model.para[model.para$name == "psiX1X1", 3],
                        model.para[model.para$name == "psiX1Xr", 3],
                        model.para[model.para$name == "psiX1X2", 3],
                        model.para[model.para$name == "psiXrXr", 3],
                        model.para[model.para$name == "psiXrX2", 3],
                        model.para[model.para$name == "psiX2X2", 3],

                        model.para[model.para$name == "betaX1Y1", 3],
                        model.para[model.para$name == "betaX1Yr", 3],
                        model.para[model.para$name == "betaX1Y2", 3],
                        model.para[model.para$name == "betaXrYr", 3],
                        model.para[model.para$name == "betaXrY2", 3],
                        model.para[model.para$name == "betaX2Y2", 3],
                        model.para[model.para$name == "betaX1M1", 3],
                        model.para[model.para$name == "betaX1Mr", 3],
                        model.para[model.para$name == "betaX1M2", 3],
                        model.para[model.para$name == "betaXrMr", 3],
                        model.para[model.para$name == "betaXrM2", 3],
                        model.para[model.para$name == "betaX2M2", 3],
                        model.para[model.para$name == "betaM1Y1", 3],
                        model.para[model.para$name == "betaM1Yr", 3],
                        model.para[model.para$name == "betaM1Y2", 3],
                        model.para[model.para$name == "betaMrYr", 3],
                        model.para[model.para$name == "betaMrY2", 3],
                        model.para[model.para$name == "betaM2Y2", 3],
                        mxSE(mediator_111, model), mxSE(mediator_rrr, model), mxSE(mediator_222, model),
                        mxSE(mediator_11r, model), mxSE(mediator_1rr, model),
                        mxSE(mediator_112, model), mxSE(mediator_1r2, model), mxSE(mediator_122, model),
                        mxSE(mediator_rr2, model), mxSE(mediator_r22, model),
                        c(mxSE(total, model)[1:3, 1], mxSE(total, model)[2:3, 2], mxSE(total, model)[3, 3]),
                        model.para[c(19, 21, 24, 20, 22, 23), 3]), 4)
    estimate_out <- data.frame(Name = paraNames, Estimate = model.est, Std.Error = model.se)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}


