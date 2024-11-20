getPosterior <- function(model, classProbs, round = NA){
  cp <- mxEval(classProbs, model)
  cp1 <- exp(cp)/sum(exp(cp))
  cp2 <- as.vector(cp1)
  cps <- diag(length(cp2))
  diag(cps) <- cp2
  subs <- model@submodels
  if (min(dim(cp)) != 1){
    stop("Class probabilities matrix must be a row or column vector.")
  }
  if (max(dim(cp)) == 1){
    stop("Class probabilities matrix must contain two or more classes.")
  }
  of <- function(num){
    return(mxEval(objective, subs[[num]]))
  }
  rl <- sapply(1:length(names(subs)), of)
  raw <- (rl %*% cps)
  tot <- 1/apply(raw, 1, sum)
  div <- matrix(rep(tot, length(cp2)), ncol=length(cp2))
  icp <- raw * div
  if (is.numeric(round)){
    icp <- round(icp, round)
  }
  return(icp)
}
