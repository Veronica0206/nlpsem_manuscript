getLCSM_JB_fixed <- function(dat, T_records, traj_var, t_var, lin_sep = 2, res_ratio = 4, starts = NA, loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA){
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  nT <- length(T_records)
  duration <- mean(dat[, paste0(t_var, nT)] - dat[, paste0(t_var, 1)])
  if (any(is.na(starts))){
    dat_traj <- dat[, paste0(traj_var, T_records)]
    dat_time <- dat[, paste0(t_var, T_records)]/duration
    eta0 <- dat_traj[, 1]
    eta1 <- (dat_traj[, nT] - dat_traj[, ceiling(nT/lin_sep)])/
      (dat_time[, nT] - dat_time[, ceiling(nT/lin_sep)])
    eta2 <- dat_traj[, 1] - dat_traj[, ceiling(nT/lin_sep)]
    growth_factor <- data.frame(eta0 = eta0, eta1 = eta1, eta2 = eta2)
    ### Initial values
    #### Mean vector
    mean0 <- c(apply(growth_factor, 2, mean))
    #### Var-cov matrix
    psi0 <- cov(growth_factor)
    ### Rate approach to asymptotic level
    slp <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      slp[j] <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[, j + 1] - dat_traj[, j])/
                                                                  (time_delta = dat_time[, j + 1] - dat_time[, j])),
                              na.action = na.exclude)$coefficients)
    }
    ratio <- rep(0, length(slp) - 1)
    for (j in 1:(length(slp) - 1)){
      ratio[j] <- slp[j + 1] - slp[j]
    }
    gamma <- mean(ratio)
    #### Decide the initial value of the residual variance
    residuals <- var(eta0)/res_ratio
    starts <- c(mean0, gamma, psi0[row(psi0) >= col(psi0)], residuals)
  }
  ### Define manifest variables
  manifests <- paste0(traj_var, T_records)
  ### Define latent variables
  latents <- c("eta0", "eta1", "eta2", paste0("ly", 1:nT), paste0("dy", 2:nT))
  outDef <- midTime <- outLag <- CHG_b <- CHG_i <- slp_est <- slp_se <- outLoads2 <- list()
  m_time <- apply(dat[, paste0(t_var, T_records)], 2, mean)
  m_mid_time <- rep(0, nT)
  for (j in 2:nT){
    m_mid_time[j] <- mean(dat[, paste0(t_var, j - 1)] + dat[, paste0(t_var, j)])/2
  }
  for (j in 1:nT){
    outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.T", j),
                            name = paste0("t", j))
    if (j == 1){
      midTime[[j]] <- mxAlgebra(0, name = paste0("mid_t", j))
      outLag[[j]] <- mxAlgebra(0, name = paste0("lag", j))
      outLoads2[[j]] <- mxAlgebra(0, name = paste0("L2", j))
      slp_est[[j]] <- mxAlgebra(0, name = paste0("instant_rate_est", j))
      slp_se[[j]] <- mxAlgebra(0, name = paste0("instant_rate_se", j))
      CHG_i[[j]] <- mxAlgebra(0, name = paste0("Change_in_Interval", j))
      CHG_b[[j]] <- mxAlgebra(0, name = paste0("Change_from_Baseline", j))
    }
    else{
      midTime[[j]] <- mxAlgebraFromString(paste0("(t", j , " +  t", j - 1, ")/2"), name = paste0("mid_t", j))
      outLag[[j]] <- mxAlgebraFromString(paste0("t", j , " -  t", j - 1), name = paste0("lag", j))
      outLoads2[[j]] <- mxAlgebraFromString(paste0("gamma_s * exp(gamma_s * mid_t", j, ")"), name = paste0("L2", j))
      slp_est[[j]] <- mxAlgebraFromString(paste0("mueta1 + gamma * mueta2 * exp(gamma * m_mid_t[", j, ",1])"), name = paste0("instant_rate_est", j))
      slp_se[[j]] <- mxAlgebraFromString(paste0("psi11 + gamma^2 * (exp(gamma * m_mid_t[", j,
                                                ",1]))^2 * psi22 + 2 * gamma * exp(gamma * m_mid_t[", j, ",1]) * psi12"), name = paste0("instant_rate_se", j))
      CHG_i[[j]] <- mxAlgebraFromString(paste0("instant_rate_est", j, " * (m_time[", j, ",1] - m_time[", j - 1, ",1])"),
                                        name = paste0("Change_in_Interval", j))
      CHG_b[[j]] <- mxAlgebraFromString(paste0("Change_from_Baseline", j - 1, " + Change_in_Interval", j),
                                        name = paste0("Change_from_Baseline", j))
    }
  }
  dat_model <- cbind(dat_traj, dat_time)
  model_mx <- mxModel("Jenss Bayley Latent Change Model", type = "RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat_model, type = "raw"),
                      #### Define factor loadings from latent variables to manifests
                      mxPath(from = "eta0", to = "ly1", arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1", to = paste0("dy", 2:nT), arrows = 1, free = F, values = 1),
                      mxPath(from = "eta2", to = paste0("dy", 2:nT), arrows = 1, free = F, values = 0,
                             labels = paste0("L2", 2:nT, "[1,1]")),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents[1:3], arrows = 1, free = T, values = starts[1:3],
                             labels = c("mueta0", "mueta1_s", "mueta2")),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents[1:3], to = latents[1:3], arrows = 2,
                             connect = "unique.pairs", free = T,
                             values = starts[5:10],
                             labels = c("psi00", "psi01_s", "psi02",
                                        "psi11_s", "psi12_s", "psi22")),
                      #### Define latent true scores
                      mxPath(from = paste0("ly", 1:nT), to = manifests, arrows = 1, free = F, values = 1),
                      #### Define path from latent instantaneous rate of change at each measurement to true scores
                      mxPath(from = paste0("dy", 2:nT), to = paste0("ly", 2:nT), arrows = 1, free = F, values = 0,
                             labels = paste0("lag", 2:nT, "[1,1]")),
                      #### Define autoregressive paths
                      mxPath(from = paste0("ly", 1:(nT - 1)), to = paste0("ly", 2:nT),
                             arrows = 1, free = F, values = 1),
                      #### Define residual variances
                      mxPath(from = manifests, to = manifests, arrows = 2, free = T, values = starts[11],
                             labels = paste0("residuals")),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = starts[4],
                               labels = paste0("gamma_s"), name = paste0("GAMMA")),
                      mxMatrix("Full", 1, 1, free = F, values = duration, name = "duration"),
                      mxAlgebra(mueta1_s/duration, name = "mueta1"),
                      mxAlgebra(gamma_s/duration, name = "gamma"),
                      mxAlgebra(psi01_s/duration, name = "psi01"),
                      mxAlgebra(psi11_s/duration^2, name = "psi11"),
                      mxAlgebra(psi12_s/duration, name = "psi12"),
                      mxAlgebra(rbind(mueta0, mueta1, mueta2), name = "mean0"),
                      mxAlgebra(rbind(cbind(psi00, psi01, psi02),
                                      cbind(psi01, psi11, psi12),
                                      cbind(psi02, psi12, psi22)), name = "psi0"),
                      outDef, midTime, outLag, outLoads2, slp_est, slp_se, CHG_i, CHG_b,

                      mxMatrix("Full", length(T_records), 1, free = F, values = m_time, name = "m_time"),
                      mxMatrix("Full", length(T_records), 1, free = F, values = m_mid_time, name = "m_mid_t"))
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(original){
    model.para <- summary(model)$parameters[, c(1, 5, 6)]
    slp_est.est <- slp_se.est <- slp_est.se <- slp_se.se <- CHG_i.est <- CHG_b.est <- slp.se <- CHG_i.se <- CHG_b.se <- rep(0, nT - 1)
    for (j in 2:nT){
      slp_est.est[j - 1] <- mxEvalByName(paste0("instant_rate_est", j), model = model)
      slp_se.est[j - 1] <- mxEvalByName(paste0("instant_rate_se", j), model = model)
      CHG_i.est[j - 1] <- mxEvalByName(paste0("Change_in_Interval", j), model = model)
      CHG_b.est[j - 1] <- mxEvalByName(paste0("Change_from_Baseline", j), model = model)
      slp_est.se[j - 1] <- mxSE(paste0("instant_rate_est", j), model = model, forceName = T)
      slp_se.se[j - 1] <- mxSE(paste0("instant_rate_se", j), model = model, forceName = T)
      CHG_i.se[j - 1] <- mxSE(paste0("Change_in_Interval", j), model = model, forceName = T)
      CHG_b.se[j - 1] <- mxSE(paste0("Change_from_Baseline", j), model = model, forceName = T)
    }
    model.est <- round(c(model$mean0$result, mxEval(gamma, model = model), model$psi0$result[row(model$psi0$result) >= col(model$psi0$result)],
                         slp_est.est, slp_se.est, CHG_i.est, CHG_b.est, model.para[1, 2]), 4)
    mean.se <- mxSE(mean0, model)
    psi.se <- mxSE(psi0, model)
    model.se <- round(c(mean.se, mxSE(gamma, model), psi.se[row(psi.se) >= col(psi.se)], slp_est.se, slp_se.se, CHG_i.se, CHG_b.se,
                        model.para[1, 3]), 4)
    estimate_out <- data.frame(Name = paraNames,
                               Estimate = model.est, SE = model.se)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}
