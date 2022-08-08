#         trace("MARSSharveyobsFI", edit=T, where = MARSSFisherI)
library(doParallel)
library(bigstatsr)
#function (MLEobj) 
MARSSharveyobsFI <- function(MLEobj) {
  paramvector <- MARSSvectorizeparam(MLEobj)
  par.names <- names(paramvector)
  num.p <- length(paramvector)
  if (num.p == 0) {
    return(obsFI = matrix(0, 0, 0))
  }
  condition.limit <- 1e+10
  condition.limit.Ft <- 1e+05
  MODELobj <- MLEobj$marss
  n <- dim(MODELobj$data)[1]
  TT <- dim(MODELobj$data)[2]
  m <- dim(MODELobj$fixed$x0)[1]
  YM <- matrix(as.numeric(!is.na(MODELobj$data)), n, TT)
  y <- MODELobj$data
  y[YM == 0] <- 0
  if (MODELobj$tinitx == 1) {
    init.state <- "x10"
  }
  else {
    init.state <- "x00"
  }
  msg <- NULL
  I.m <- diag(1, m)
  I.n <- diag(1, n)
  kf <- try(MARSSkfss(MLEobj), silent = TRUE)
  if (inherits(kf, "try-error")) 
    stop("Stopped in MARSSharveyobsFI(). MARSSkfss does not run for this model. Try a numerical Hessian?", 
         call. = FALSE)
  xtt <- kf$xtt
  xtt1 <- kf$xtt1
  Vtt <- kf$Vtt
  Vtt1 <- kf$Vtt1
  vt <- kf$Innov
  Ft <- kf$Sigma
  dxtt1 <- matrix(0, m, 1)
  dVtt1 <- matrix(0, m, m)
  dxt1t1 <- matrix(0, m, num.p)
  dVt1t1 <- array(0, dim = c(m, m, num.p))
  dvt <- matrix(0, n, num.p)
  dFt <- array(0, dim = c(n, n, num.p))
  obsFI <- matrix(0, num.p, num.p)
  rownames(obsFI) <- par.names
  colnames(obsFI) <- par.names
  model.elem <- attr(MODELobj, "par.names")
  dims <- attr(MODELobj, "model.dims")
  time.varying <- c()
  for (el in model.elem) {
    if ((dim(MODELobj$free[[el]])[3] != 1) || (dim(MODELobj$fixed[[el]])[3] != 
                                               1)) {
      time.varying <- c(time.varying, el)
    }
  }
  pari <- parmat(MLEobj, t = 1)
  Z <- pari$Z
  A <- pari$A
  B <- pari$B
  U <- pari$U
  x0 <- pari$x0
  R <- tcrossprod(pari$H %*% pari$R, pari$H)
  Q <- tcrossprod(pari$G %*% pari$Q, pari$G)
  V0 <- tcrossprod(pari$L %*% pari$V0, pari$L)
  dpari <- dparmat(MLEobj, t = 1)
  time.varying.par <- time.varying[time.varying %in% names(dpari)]
  dpar0 <- list()
  for (el in model.elem) dpar0[[el]] <- matrix(0, dims[[el]][1], 
                                               dims[[el]][2])
  for (t in 1:TT) {
    print(c("TIME",t))
    if (length(time.varying) != 0) {
      pari[time.varying] <- parmat(MLEobj, time.varying, 
                                   t = t)
      Z <- pari$Z
      A <- pari$A
      B <- pari$B
      U <- pari$U
      dpari[time.varying.par] <- dparmat(MLEobj, time.varying.par, 
                                         t = t)
    }
    pcntr <- 0
    for (el in model.elem) {
      dp <- dpar0
      p <- length(MLEobj$par[[el]])
      if (p == 0) 
        next
      for (ip in 1:p) {
        pcntr <- pcntr + 1
        dp[[el]] <- array(dpari[[el]][, , ip], dim = dims[[el]][1:2])
        dHRH <- tcrossprod(dp$H %*% pari$R, pari$H) + 
          tcrossprod(pari$H %*% dp$R, pari$H) + tcrossprod(pari$H %*% 
                                                             pari$R, dp$H)
        dGQG <- tcrossprod(dp$G %*% pari$Q, pari$G) + 
          tcrossprod(pari$G %*% dp$Q, pari$G) + tcrossprod(pari$G %*% 
                                                             pari$Q, dp$G)
        dLV0L <- tcrossprod(dp$L %*% pari$V0, pari$L) + 
          tcrossprod(pari$L %*% dp$V0, pari$L) + tcrossprod(pari$L %*% 
                                                              pari$V0, dp$L)
        if (any(YM[, t] == 0)) {
          Mt <- I.n
          Mt[YM[, t] == 0, ] <- 0
          I.2 <- I.n - Mt
          Zt <- Mt %*% Z
          dp$Zt <- Mt %*% dp$Z
          dp$At <- Mt %*% dp$A
          dHRHt <- Mt %*% dHRH %*% Mt
        }
        else {
          Zt <- Z
          dp$Zt <- dp$Z
          dHRHt <- dHRH
          dp$At <- dp$A
        }
        if (t == 1) {
          if (init.state == "x00") {
            dxtt1 <- dp$B %*% x0 + B %*% dp$x0 + dp$U
            dVtt1 <- tcrossprod(dp$B %*% V0, B) + tcrossprod(B %*% 
                                                               dp$V0, B) + tcrossprod(B %*% V0, dp$B) + 
              dGQG
          }
          if (init.state == "x10") {
            dxtt1 <- dp$x0
            dVtt1 <- dp$V0
          }
        }
        else {
          dxtt1 <- dp$B %*% xtt[, t - 1, drop = FALSE] + 
            B %*% dxt1t1[, pcntr, drop = FALSE] + dp$U
          dVtt1 <- tcrossprod(dp$B %*% Vtt[, , t - 1], 
                              B) + tcrossprod(B %*% dVt1t1[, , pcntr], 
                                              B) + tcrossprod(B %*% Vtt[, , t - 1], dp$B) + 
            dGQG
        }
        if (m != 1) 
          dVtt1 <- symm(dVtt1)
        dvt[, pcntr] <- -Zt %*% dxtt1 - dp$Zt %*% xtt1[, 
                                                       t, drop = FALSE] - dp$At
        dFt[, , pcntr] <- tcrossprod(dp$Zt %*% Vtt1[, 
                                                    , t], Zt) + tcrossprod(Zt %*% dVtt1, Zt) + 
          tcrossprod(Zt %*% Vtt1[, , t], dp$Zt) + dHRHt
        if (n == 1) {
          Ftinv <- pcholinv(matrix(Ft[, , t], 1, 1))
        }
        else {
          Ftinv <- pcholinv(Ft[, , t])
          Ftinv <- symm(Ftinv)
        }
        dxt1t1[, pcntr] <- dxtt1 + tcrossprod(dVtt1, 
                                              Zt) %*% Ftinv %*% vt[, t, drop = FALSE] + tcrossprod(Vtt1[, 
                                                                                                        , t], dp$Zt) %*% Ftinv %*% vt[, t, drop = FALSE] - 
          tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% dFt[, 
                                                        , pcntr] %*% Ftinv %*% vt[, t, drop = FALSE] + 
          tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% dvt[, 
                                                        pcntr, drop = FALSE]
        dVt1t1[, , pcntr] <- dVtt1 - tcrossprod(dVtt1, 
                                                Zt) %*% Ftinv %*% Zt %*% Vtt1[, , t] - tcrossprod(Vtt1[, 
                                                                                                       , t], dp$Zt) %*% Ftinv %*% Zt %*% Vtt1[, , 
                                                                                                                                              t] + tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% 
          dFt[, , pcntr] %*% Ftinv %*% Zt %*% Vtt1[, 
                                                   , t] - tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% 
          dp$Zt %*% Vtt1[, , t] - tcrossprod(Vtt1[, , 
                                                  t], Zt) %*% Ftinv %*% Zt %*% dVtt1
        if (m != 1) 
          dVt1t1[, , pcntr] <- symm(dVt1t1[, , pcntr])
      }
    }
    mat3<-FBM(num.p,num.p)
    print(c("PARA BEGIN NEW"))
    print(Sys.time())
    registerDoParallel(25)
    tmp3<-foreach(i=1:num.p,.combine="c") %:%
      foreach(j=1:num.p,.combine = "c") %dopar% {
        mat3[i,j] <- sum(diag(Ftinv %*% dFt[, , i] %*% 
                                Ftinv %*% dFt[, , j]))/2 + t(dvt[, i, drop = FALSE]) %*% 
          Ftinv %*% dvt[, j, drop = FALSE]
        NULL
      }
    stopImplicitCluster()
    tmp<-mat3[]
    print(c("PARA FINISHED"))
    print(Sys.time())
    obsFI <- obsFI + tmp
  }
  return(obsFI)
}
