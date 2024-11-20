getBLSGM_Fixed <- function(dat, T_records, traj_var, t_var, res_ratio = 4, starts = NA, loc = 1, scale = 0.25, extraTries = NA, original = T,
                           paraNames = NA){
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  nT <- length(T_records)
  if (any(is.na(starts))){
    dat_traj <- dat[, paste0(traj_var, T_records)]
    dat_time <- dat[, paste0(t_var, T_records)]
    slp <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      slp[j] <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[, j + 1] - dat_traj[, j])/
                                                                  (time_delta = dat_time[, j + 1] - dat_time[, j])),
                              na.action = na.exclude)$coefficients)
    }
    slp_diff <- diff(slp)
    for (ratio in seq(2, 0, -0.05)){
      knot_candidate <- which(abs(slp_diff) > ratio * abs(mean((dat_traj[, nT] - dat_traj[, 1])/(dat_time[, nT] - dat_time[, 1]), na.rm = T)))
      if (length(knot_candidate) >= 2){
        #### Decide the initial value of the knot mean
        gamma <- apply(dat_time[, knot_candidate + 1], 1, mean, na.rm = T)
        break
      }
    }
    for (j in 1:(nT - 1)){
      if (I(mean(gamma) > mean(dat_time[, j]) & (mean(gamma) < mean(dat_time[, j + 1])))){
        gamma_f <- j
        gamma_c <- j + 1
      }
    }
    eta0 <- dat_traj[, 1]
    eta1 <- (dat_traj[, gamma_f] - dat_traj[, 1])/(dat_time[, gamma_f] - dat_time[, 1])
    eta2 <- (dat_traj[, nT] - dat_traj[, gamma_c])/(dat_time[, nT] - dat_time[, gamma_c])
    growth_factor <- data.frame(eta0 = eta0, eta1 = eta1, eta2 = eta2, gamma = gamma)

    ### Initial values in the original parameter space
    #### Mean vector
    mean0 <- apply(growth_factor, 2, mean)
    #### Var-cov matrix
    psi0 <- cov(growth_factor[, 1:3])
    #### Decide the initial value of the residual variance
    residuals <- var(eta0)/res_ratio
    ### Transformed matrices obtained by multivariate Delta method
    #### For mean vector
    func0 <- matrix(c(1, mean0[4], 0,
                      0, 0.5, 0.5,
                      0, -0.5, 0.5), nrow = 3, byrow = T)
    #### For var-cov matrix
    grad0 <- matrix(c(1, mean0[4], 0,
                      0, 0.5, 0.5,
                      0, -0.5, 0.5), nrow = 3, byrow = T)

    mean0.s <- func0 %*% mean0[1:3]
    psi0.s <- grad0 %*% psi0 %*% t(grad0)
    starts <- c(mean0.s[1:3], mean0[4], psi0.s[row(psi0.s) >= col(psi0.s)], residuals)
  }
  ### Define manifest variables
  manifests <- paste0(traj_var, T_records)
  ### Define latent variables
  latents <- c("eta0s", "eta1s", "eta2s")
  outDef <- list(); outLoads1 <- list(); outLoads2 <- list()
  for(j in T_records){
    outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var, j), name = paste0("t", j))
    outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " - mug"), name = paste0("L1", j))
    outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " - mug)"), name = paste0("L2", j))
  }
  ### Create a mxModel object
  model_mx <- mxModel("Bilinear Spline Growth Model with a Fixed Knot", type = "RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define factor loading from latent variables to manifests
                      mxPath(from = "eta0s", to = manifests, arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1s", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L1", T_records, "[1,1]")),
                      mxPath(from = "eta2s", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L2", T_records, "[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = manifests, to = manifests, arrows = 2, free = T, values = starts[11],
                             labels = "residuals"),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents, arrows = 1, free = T, values = starts[1:3],
                             labels = c("mueta0s", "mueta1s", "mueta2s")),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[5:10],
                             labels = c("psi0s0s", "psi0s1s", "psi0s2s",
                                        "psi1s1s", "psi1s2s", "psi2s2s")),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = starts[4],
                               labels = "muknot", name = "mug"),
                      outDef, outLoads1, outLoads2,
                      mxAlgebra(rbind(mueta0s, mueta1s, mueta2s), name = "mean_s"),
                      mxAlgebra(rbind(cbind(psi0s0s, psi0s1s, psi0s2s),
                                      cbind(psi0s1s, psi1s1s, psi1s2s),
                                      cbind(psi0s2s, psi1s2s, psi2s2s)), name = "psi_s"),
                      mxAlgebra(rbind(cbind(1, -mug, mug),
                                      cbind(0, 1, -1),
                                      cbind(0, 1, 1)), name = "func"),
                      mxAlgebra(rbind(cbind(1, -mug, mug),
                                      cbind(0, 1, -1),
                                      cbind(0, 1, 1)), name = "grad"),
                      mxAlgebra(func %*% mean_s, name = "mean"),
                      mxAlgebra(grad %*% psi_s %*% t(grad), name = "psi"))
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(original){
    model.para <- summary(model)$parameters[, c(1, 5, 6)]
    model.est <- round(c(model$mean$result, model.para[model.para$name == "muknot", 2],
                         model$psi$result[row(model$psi$result) >= col(model$psi$result)],
                         model.para[model.para$name == "residuals", 2]), 4)
    mean.se <- mxSE(mean, model)
    psi.se <- mxSE(psi, model)
    model.se <- round(c(mean.se, model.para[model.para$name == "muknot", 3],
                        psi.se[row(psi.se) >= col(psi.se)],
                        model.para[model.para$name == "residuals", 3]), 4)
    estimate_out <- data.frame(Name = paraNames, Estimate = model.est, SE = model.se)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}


