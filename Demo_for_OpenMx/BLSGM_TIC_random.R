getBLSGM_TIC_Random <- function(dat, T_records, traj_var, t_var, growth_cov, res_ratio = 4, gv_adjust = 1, rho_gamma = rep(0.3, 3), starts = NA,
                                loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA){
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  if (any(is.na(starts))){
    dat_traj <- dat[, paste0(traj_var, T_records)]
    dat_time <- dat[, paste0(t_var, T_records)]
    dat_covariate <- dat[, growth_cov]
    ### Obtain estimates from the reduced model
    reduced_model <- getBLSGM_TIC_Fixed(dat = dat, T_records = T_records, traj_var = traj_var, t_var = t_var, growth_cov = growth_cov,
                                        res_ratio = res_ratio, original = F)
    ### Obtain knot variance
    mean_time <- apply(dat_time, 2, mean)
    mug <- reduced_model$mug$values
    ### Obtain knot variance and path coefficients to knot
    mean_time <- apply(dat_time, 2, mean)
    mug <- reduced_model$mug$values
    for (j in 1:length(mean_time)){
      if (mean_time[j] <= mug & mean_time[j + 1] >= mug){
        BetaG <- as.numeric(lm(knot ~ ., data = data.frame(cbind(knot = apply(dat_time[, c(j, j + 1)], 1, mean, na.rm = T) * gv_adjust, dat_covariate)),
                               na.action = na.exclude)$coefficients[-1])
        psigg <- var(apply(dat_time[, c(j, j + 1)], 1, mean, na.rm = T)) * gv_adjust -
          t(BetaG) %*% diag(apply(dat_covariate, 2, var, na.rm = T)) %*% BetaG
        stop
      }
    }

    ### Initial values in the original parameter space
    mean0 <- c(reduced_model$mean$result, mug)
    reduced_psi0 <- reduced_model$psi$result
    psi0 <- matrix(c(reduced_psi0[1, ], rho_gamma[1] * sqrt(reduced_psi0[1, 1] * psigg),
                     reduced_psi0[2, ], rho_gamma[2] * sqrt(reduced_psi0[2, 2] * psigg),
                     reduced_psi0[3, ], rho_gamma[3] * sqrt(reduced_psi0[3, 3] * psigg),
                     rho_gamma[1] * sqrt(reduced_psi0[1, 1] * psigg), rho_gamma[2] * sqrt(reduced_psi0[2, 2] * psigg),
                     rho_gamma[3] * sqrt(reduced_psi0[3, 3] * psigg), psigg), nrow = 4)
    beta0 <- matrix(c(reduced_model$beta$result[1, ], reduced_model$beta$result[2, ],
                      reduced_model$beta$result[3, ], BetaG), nrow = 4, byrow = T)
    ### Transformed matrices obtained by multivariate Delta method
    #### For mean vector
    func0 <- matrix(c(1, mean0[4], 0, 0,
                      0, 0.5, 0.5, 0,
                      0, -0.5, 0.5, 0,
                      0, 0, 0, 1), nrow = 4, byrow = T)
    #### For var-cov matrix
    grad0 <- matrix(c(1, mean0[4], 0, mean0[2],
                      0, 0.5, 0.5, 0,
                      0, -0.5, 0.5, 0,
                      0, 0, 0, 1), nrow = 4, byrow = T)

    mean0.s <- func0[1:3, 1:3] %*% mean0[1:3]
    psi0.s <- grad0 %*% psi0 %*% t(grad0)
    beta.s <- grad0 %*% beta0
    starts.y <- c(mean0.s[1:3], mean0[4], psi0.s[row(psi0.s) >= col(psi0.s)], reduced_model$S$values[1, 1])
    starts.x <- c(apply(dat_covariate, 2, mean, na.rm = T),
                  var(dat_covariate)[row(var(dat_covariate)) >= col(var(dat_covariate))])
    starts.beta <- c(beta.s[1, ], beta.s[2, ], beta.s[3, ], beta.s[4, ])
    starts <- list(starts.y, starts.x, starts.beta)
  }
  ### Define manifest variables
  manifests <- paste0(traj_var, T_records)
  ### Define latent variables
  latents <- c("eta0s", "eta1s", "eta2s", "delta")
  outDef <- list(); outLoads1 <- list(); outLoads2 <- list(); outLoads3 <- list()
  for(j in T_records){
    outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var, j), name = paste0("t", j))
    outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " - mug"), name = paste0("L1", j))
    outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " - mug)"), name = paste0("L2", j))
    outLoads3[[j]] <- mxAlgebraFromString(paste0("-mueta2s * (t", j, " - mug)/abs(t", j,
                                                 " - mug) - mueta2s"), name = paste0("L3", j))
  }
  ### Create a mxModel object
  model_mx <- mxModel("Bilinear Spline Growth Model with a Random Knot and Growth Covariates", type = "RAM",
                      manifestVars = c(manifests, growth_cov), latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define factor loadings from latent variables to manifests
                      mxPath(from = "eta0s", to = manifests, arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1s", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L1", T_records, "[1,1]")),
                      mxPath(from = "eta2s", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L2", T_records, "[1,1]")),
                      mxPath(from = "delta", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L3", T_records, "[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = manifests, to = manifests, arrows = 2, free = T, values = starts[[1]][15],
                             labels = "residuals"),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents[1:3], arrows = 1, free = T, values = starts[[1]][1:3],
                             labels = c("mueta0s", "mueta1s", "mueta2s")),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents, to = latents, arrows = 2,
                             connect = "unique.pairs", free = T,
                             values = starts[[1]][5:14],
                             labels = c("psi0s0s", "psi0s1s", "psi0s2s", "psi0sg", "psi1s1s",
                                        "psi1s2s", "psi1sg", "psi2s2s", "psi2sg", "psigg")),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = starts[[1]][4],
                               labels = "muknot", name = "mug"),
                      #### Include time-invariant covariates
                      ##### Means
                      mxPath(from = "one", to = growth_cov, arrows = 1, free = T,
                             values = starts[[2]][1:length(growth_cov)], labels = paste0("mux", 1:length(growth_cov))),
                      mxPath(from = growth_cov, to = growth_cov, connect = "unique.pairs",
                             arrows = 2, free = T,
                             values = starts[[2]][-1:-length(growth_cov)],
                             labels = paste0("phi", 1:(length(growth_cov) * (length(growth_cov) + 1)/2))),
                      ##### Regression coefficients
                      mxPath(from = growth_cov, to = latents[1], arrows = 1, free = T,
                             values = starts[[3]][1:length(growth_cov)],
                             labels = paste0("beta0", 1:length(growth_cov))),
                      mxPath(from = growth_cov, to = latents[2], arrows = 1, free = T,
                             values = starts[[3]][(2 + 1):(2 + length(growth_cov))],
                             labels = paste0("beta1", 1:length(growth_cov))),
                      mxPath(from = growth_cov, to = latents[3], arrows = 1, free = T,
                             values = starts[[3]][(2 * 2 + 1):(2 * 2 + length(growth_cov))],
                             labels = paste0("beta2", 1:length(growth_cov))),
                      mxPath(from = growth_cov, to = latents[4], arrows = 1, free = T,
                             values = starts[[3]][(3 * 2 + 1):(3 * 2 + length(growth_cov))],
                             labels = paste0("betaG", 1:length(growth_cov))),
                      outDef, outLoads1, outLoads2, outLoads3,
                      mxAlgebra(rbind(mueta0s, mueta1s, mueta2s, mug), name = "mean_s"),
                      mxAlgebra(rbind(cbind(psi0s0s, psi0s1s, psi0s2s, psi0sg),
                                      cbind(psi0s1s, psi1s1s, psi1s2s, psi1sg),
                                      cbind(psi0s2s, psi1s2s, psi2s2s, psi2sg),
                                      cbind(psi0sg, psi1sg, psi2sg, psigg)), name = "psi_s"),
                      mxMatrix("Full", 4, length(growth_cov), free = T, values = starts[[3]],
                               labels = c(paste0("beta0", 1:length(growth_cov)),
                                          paste0("beta1", 1:length(growth_cov)),
                                          paste0("beta2", 1:length(growth_cov)),
                                          paste0("betaG", 1:length(growth_cov))), byrow = T, name = "beta_s"),
                      mxMatrix("Full", 1, length(growth_cov), free = T, values = starts[[2]][1:length(growth_cov)],
                               labels = paste0("mux", 1:length(growth_cov)), byrow = T, name = "BL_mean"),
                      mxAlgebra(rbind(cbind(1, -mug, mug, 0),
                                      cbind(0, 1, -1, 0),
                                      cbind(0, 1, 1, 0),
                                      cbind(0, 0, 0, 1)), name = "func"),
                      mxAlgebra(rbind(cbind(1, -mug, mug, 0),
                                      cbind(0, 1, -1, 0),
                                      cbind(0, 1, 1, 0),
                                      cbind(0, 0, 0, 1)), name = "grad"),
                      mxAlgebra(func %*% (mean_s + beta_s %*% t(BL_mean)), name = "mean"),
                      mxAlgebra(grad %*% psi_s %*% t(grad), name = "psi"),

                      mxAlgebra(grad %*% beta_s, name = "beta"))
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(original){
    model.para <- summary(model)$parameters[, c(1, 5, 6)]
    model.est <- round(c(model$mean$result,
                         model$psi$result[row(model$psi$result) >= col(model$psi$result)],
                         c(model$beta$result),
                         model.para[grep("mux", model.para$name), 2],
                         model.para[grep("phi", model.para$name), 2],
                         model.para[model.para$name == "residuals", 2]), 4)
    mean.se <- mxSE(mean, model)

    psi.se <- mxSE(psi, model)
    beta.se <- mxSE(beta, model)
    model.se <- round(c(mean.se,
                        psi.se[row(psi.se) >= col(psi.se)], c(beta.se),
                        model.para[grep("mux", model.para$name), 3],
                        model.para[grep("phi", model.para$name), 3],
                        model.para[model.para$name == "residuals", 3]), 4)
    estimate_out <- data.frame(Name = paraNames, Estimate = model.est, SE = model.se)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}
