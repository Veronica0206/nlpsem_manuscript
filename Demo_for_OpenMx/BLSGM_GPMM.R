getGPMM_BLSGM <- function(dat, T_records, nClass, traj_var, t_var, growth_cov, res_ratio = rep(4, 2), starts = NA, prop_starts = rep(0.5, 2),
                          loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA){
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  else{
    if (any(is.na(starts))){
      starts <- list()
      FMM <- getFMM_BLSGM(dat = dat, T_records = T_records, nClass = nClass, traj_var = traj_var, t_var = t_var, res_ratio = res_ratio,
                          prop_starts = prop_starts, original = F)
      model.para <- summary(FMM)$parameters[, c(1, 5, 6)]
      FMM_post <- getPosterior(model = FMM, classProbs = FMM$weights$values, round = 4)
      dat$label <- apply(FMM_post, 1, which.max)
      for (k in 1:nClass){
        subdat <- dat[dat$label == k, ]
        submodel <- getBLSGM_TIC_Fixed(dat = subdat, T_records = T_records, traj_var = traj_var, t_var = t_var, growth_cov = growth_cov,
                                       res_ratio = res_ratio[k], original = F)
        starts[[length(starts) + 1]] <- list(c(submodel$mean_s$result, submodel$mug$values,
                                               submodel$psi_s$result[row(submodel$psi_s$result) >= col(submodel$psi_s$result)],
                                               submodel$S$values[1, 1]),
                                             c(submodel@output$estimate[grep("mux", names(submodel@output$estimate))],
                                               submodel@output$estimate[grep("phi", names(submodel@output$estimate))]),
                                             c(submodel$beta_s$values[1, ], submodel$beta_s$values[2, ], submodel$beta_s$values[3, ]))
      }
    }
    prop_update <- table(dat$label)/sum(table(dat$label))
    w_starts <- c(log(prop_update/prop_update[1] * exp(1)))
    names(w_starts) <- NULL
  }
  ### Define manifest variables
  manifests <- paste0(traj_var, T_records)
  ### Define latent variables
  latents <- c("eta0s", "eta1s", "eta2s")
  outDef <- list(); outLoads1 <- list(); outLoads2 <- list()
  class.list <- list()
  for (k in 1:nClass){
    for(j in T_records){
      outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.T", j),
                              name = paste0("t", j))
      outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " -c", k, "mug"), name = paste0("c", k, "L1", j))
      outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " -c", k, "mug)"),
                                            name = paste0("c", k, "L2", j))
    }
    ### Create a mxModel object
    class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                               manifestVars = c(manifests, growth_cov), latentVars = latents,
                               #### Define factor loadings from latent variables to manifests
                               mxPath(from = "eta0s", to = manifests, arrows = 1, free = F, values = 1),
                               mxPath(from = "eta1s", to = manifests, arrows = 1, free = F, values = 0,
                                      labels = paste0("c", k, "L1", T_records, "[1,1]")),
                               mxPath(from = "eta2s", to = manifests, arrows = 1, free = F, values = 0,
                                      labels = paste0("c", k, "L2", T_records, "[1,1]")),
                               #### Define the variances of residuals
                               mxPath(from = manifests, to = manifests, arrows = 2, free = T, values = starts[[k]][[1]][11],
                                      labels = paste0("c", k, "residuals")),
                               #### Define means of latent variables
                               mxPath(from = "one", to = latents, arrows = 1, free = T, values = starts[[k]][[1]][1:3],
                                      labels = paste0("c", k, c("mueta0s", "mueta1s", "mueta2s"))),
                               #### Define var-cov matrix of latent variables
                               mxPath(from = latents, to = latents, arrows = 2,
                                      connect = "unique.pairs", free = T,
                                      values = starts[[k]][[1]][c(5:10)],
                                      labels = paste0("c", k, c("psi0s0s", "psi0s1s", "psi0s2s",
                                                                "psi1s1s", "psi1s2s", "psi2s2s"))),
                               #### Add additional parameter and constraints
                               mxMatrix("Full", 1, 1, free = T, values = starts[[k]][[1]][4],
                                        labels = paste0("c", k, "muknot"), name = paste0("c", k, "mug")),
                               #### Include time-invariant covariates
                               ##### Means
                               mxPath(from = "one", to = growth_cov, arrows = 1, free = T,
                                      values = starts[[k]][[2]][1:length(growth_cov)],
                                      labels = paste0("c", k, "mux", 1:length(growth_cov))),
                               mxPath(from = growth_cov, to = growth_cov, connect = "unique.pairs",
                                      arrows = 2, free = T,
                                      values = starts[[k]][[2]][-1:-length(growth_cov)],
                                      labels = paste0("c", k, "phi", 1:(length(growth_cov) * (length(growth_cov) + 1)/2))),
                               ##### Regression coefficients
                               mxPath(from = growth_cov, to = latents[1], arrows = 1, free = T,
                                      values = starts[[k]][[3]][1:length(growth_cov)],
                                      labels = paste0("c", k, "beta0", 1:length(growth_cov))),
                               mxPath(from = growth_cov, to = latents[2], arrows = 1, free = T,
                                      values = starts[[k]][[3]][(2 + 1):(2 + length(growth_cov))],
                                      labels = paste0("c", k, "beta1", 1:length(growth_cov))),
                               mxPath(from = growth_cov, to = latents[3], arrows = 1, free = T,
                                      values = starts[[k]][[3]][(2 * 2 + 1):(2 * 2 + length(growth_cov))],
                                      labels = paste0("c", k, "beta2", 1:length(growth_cov))),
                               outDef, outLoads1, outLoads2,
                               mxAlgebraFromString(paste0("rbind(c", k, "mueta0s, c", k, "mueta1s, c", k, "mueta2s)"),
                                                   name = paste0("c", k, "mean_s")),
                               mxAlgebraFromString(paste0("rbind(cbind(c", k, "psi0s0s, c", k, "psi0s1s, c", k, "psi0s2s),",
                                                          "cbind(c", k, "psi0s1s, c", k, "psi1s1s, c", k, "psi1s2s),",
                                                          "cbind(c", k, "psi0s2s, c", k, "psi1s2s, c", k, "psi2s2s))"),
                                                   name = paste0("c", k, "psi_s")),
                               mxMatrix("Full", 3, length(growth_cov), free = T, values = starts[[k]][[3]],
                                        labels = c(paste0("c", k, "beta0", 1:length(growth_cov)),
                                                   paste0("c", k, "beta1", 1:length(growth_cov)),
                                                   paste0("c", k, "beta2", 1:length(growth_cov))), byrow = T, name = paste0("c", k, "beta_s")),
                               mxMatrix("Full", 1, length(growth_cov), free = T, values = starts[[k]][[2]][1:length(growth_cov)],
                                        labels = paste0("c", k, "mux", 1:length(growth_cov)), byrow = T, name = paste0("c", k, "BL_mean")),
                               mxAlgebraFromString(paste0("rbind(cbind(", "1,", "-c", k, "muknot,", "c", k, "muknot),",
                                                          "cbind(0, 1, -1), cbind(0, 1, 1))"),
                                                   name = paste0("c", k, "func")),
                               mxAlgebraFromString(paste0("rbind(cbind(", "1,", "-c", k, "muknot,", "c", k, "muknot),",
                                                          "cbind(0, 1, -1), cbind(0, 1, 1))"),
                                                   name = paste0("c", k, "grad")),
                               mxAlgebraFromString(paste0("c", k, "func %*% (c", k, "mean_s + c", k, "beta_s %*% t(c", k, "BL_mean))"), name = paste0("c", k, "mean")),
                               mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "psi_s %*% t(c", k, "grad)"),
                                                   name = paste0("c", k, "psi")),
                               mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "beta_s"),
                                                   name = paste0("c", k, "beta")),
                               mxFitFunctionML(vector = T))
  }
  ### Make the class proportion matrix, fixing one parameter at a non-zero constant (one)
  classP <- mxMatrix("Full", nClass, 1, free = c(F, rep(T, nClass - 1)), values = w_starts,
                     labels = paste0("w", 1:nClass), name = "weights")
  algebraObjective <- mxExpectationMixture(paste0("Class", 1:nClass), weights = "weights", scale = "softmax")
  objective <- mxFitFunctionML()
  model_mx <- mxModel("Growth Predictor Mixture Model, Bilinear Spline Growth Model with a Fixed Knot",
                      mxData(observed = dat, type = "raw"), class.list, classP, algebraObjective, objective)
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if (original){
    model.para <- summary(model)$parameters[, c(1, 5, 6)]
    ### Estimates of parameters in each latent class
    model.est <- model.se <- est <- list()
    for (k in 1:nClass){
      model.est[[k]] <- round(c(mxEvalByName(paste0("c", k, "mean"), model = model@submodels[[k]]),
                                model.para[model.para$name == paste0("c", k, "muknot"), 2],
                                mxEvalByName(paste0("c", k, "psi"), model = model@submodels[[k]])[
                                  row(mxEvalByName(paste0("c", k, "psi"), model = model@submodels[[k]])) >=
                                    col(mxEvalByName(paste0("c", k, "psi"), model = model@submodels[[k]]))],
                                mxEvalByName(paste0("c", k, "beta"), model = model@submodels[[k]]),
                                model.para[grep(paste0("c", k, "mux"), model.para$name), 2],
                                model.para[grep(paste0("c", k, "phi"), model.para$name), 2],
                                model.para[model.para$name == paste0("c", k, "residuals"), 2]), 4)
      model.se[[k]] <- round(c(mxSE(paste0("Class", k, ".c", k, "mean"), model, forceName = T),
                               model.para[model.para$name == paste0("c", k, "muknot"), 3],
                               mxSE(paste0("Class", k, ".c", k, "psi"), model, forceName = T)[
                                 row(mxSE(paste0("Class", k, ".c", k, "psi"), model, forceName = T)) >=
                                   col(mxSE(paste0("Class", k, ".c", k, "psi"), model, forceName = T))],
                               mxSE(paste0("Class", k, ".c", k, "beta"), model, forceName = T),
                               model.para[grep(paste0("c", k, "mux"), model.para$name), 3],
                               model.para[grep(paste0("c", k, "phi"), model.para$name), 3],
                               model.para[model.para$name == paste0("c", k, "residuals"), 3]), 4)
      est[[k]] <- data.frame(Name = paste0("c", k, paraNames),
                             Estimate = model.est[[k]], SE = model.se[[k]])
    }
    est.weights <- data.frame(Name = paste0("p", 2:nClass),
                              Estimate = round(model.para[grep("w", model.para$name), 2], 4),
                              SE = round(model.para[grep("w", model.para$name), 3], 4))
    estimate_out <- rbind(do.call(rbind.data.frame, est), est.weights)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}
