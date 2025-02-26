sar_loglik <- function(f_sar,
                       pars,
                       area,
                       obs_richness) {
  exp_richness <- f_sar(area, pars)
  loglik <- sum(dpois(x = obs_richness, lambda = exp_richness, log = TRUE))
  return(loglik)
}

f_asymp <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- d - c * z^area)
}

f_betap <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  f <- pars[3]
  z <- pars[4]
  return(S <- d * (1 - (1 + (area/c)^z)^(-f)))
}

f_chapman <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- d * (1 - exp(-z * area)^c))
}

f_loga <- function(area, pars) {
  c <- pars[1]
  z <- pars[2]
  return(S <- c + z * log(area))
}

f_epm1 <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- c * area^(z * area^(-d)))
}

f_epm2 <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- c * area^(z - d/area))
}

f_gompertz <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- d * exp(-exp(-z * (area - c))))
}

f_koba <- function(area, pars) {
  c <- pars[1]
  z <- pars[2]
  return(S <- c * log(1 + area/z))
}

f_linear <- function(area, pars) {
  c <- pars[1]
  z <- pars[2]
  return(S <- c + z * area)
}

f_heleg <- function(area, pars) {
  c <- pars[1]
  f <- pars[2]
  z <- pars[3]
  return(S <- c/(f + area^(-z)))
}

f_monod <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  return(S <- d/(1 + c * area^(-1)))
}

f_mmf <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- d/(1 + c * area^(-z)))
}

f_negexpo <- function(area, pars) {
  d <- pars[1]
  z <- pars[2]
  return(S <- d * (1 - exp(-z * area)))
}

f_p1 <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- c * area^z * exp(-d * area))
}

f_p2 <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- c * area^z * exp(-d / area))
}

f_power <- function(area, pars) {
  c <- pars[1]
  z <- pars[2]
  return(S <- c * area^z)
}

f_powerR <- function(area, pars) {
  c <- pars[1]
  f <- pars[2]
  z <- pars[3]
  return(S <- f + c * area^z)
}

f_ratio <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- (c + z * area)/(1 + d * area))
}

f_weibull3 <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  z <- pars[3]
  return(S <- d * (1 - exp(-c * area^z)))
}

f_weibull4 <- function(area, pars) {
  c <- pars[1]
  d <- pars[2]
  f <- pars[3]
  z <- pars[4]
  return(S <- d * (1 - exp(-c * area^z))^f)
}

f_ContOne <- function(area, pars) {
  c_1 <- pars[1]
  T <- pars[2]
  z_1 <- pars[3]
  z_2 <- pars[4]
  return(S <- c_1 + (log(area) <= T) * z_1 * log(area) +
           (log(area) > T) * (z_1 * T + z_2 * (log(area) - T)))
}

f_ContTwo <- function(area, pars) {
  c_1 <- pars[1]
  T_1 <- pars[2]
  T_2 <- pars[3]
  z_1 <- pars[4]
  z_2 <- pars[5]
  z_3 <- pars[6]
  return(S <- c_1 + (log(area) <= T_1) * z_1 * log(area) +
           (log(area) > T_1) * (log(area) <= T_2) * (z_1 * T_1 + z_2 * (log(area) - T_1)) +
           (log(area) > T_2) * (z_2 * (T_2 - T_1) + z_3 * (log(area) - T_2)))
}

sar_ML <- function(f_sar,
                   initparsopt,
                   parsfix,
                   idparsopt,
                   idparsfix,
                   area,
                   obs_richness,
                   optimmethod = 'simplex',
                   num_cycles = 1,
                   tol = c(1E-6, 1E-6, 1E-6),
                   maxiter = 10000 * round((1.25)^length(idparsopt)),
                   trial_settings = c(n = 10, sd = 0.1)) {
  optimpars <- c(tol, maxiter)
  ML <- -Inf
  for(i in 0:trial_settings[1]) {
    if(i > 0) {
      initparsopt <- initparsopt * exp(rnorm(n = length(initparsopt), mean = 0,sd = trial_settings[2]))
      if(length(parsfix) > 0) {
        parsfix <- parsfix * exp(rnorm(n = length(parsfix), mean = 0,sd = trial_settings[2]))
      }
    }
    trparsopt <- initparsopt/(1 + initparsopt)
    trparsfix <- parsfix/(1 + parsfix)
    trparsfix[parsfix == Inf] = 1
    out2 <- DDD::optimizer(optimmethod = optimmethod,
                          optimpars = optimpars,
                          fun = sar_loglik_choosepar,
                          trparsopt = trparsopt,
                          trparsfix = trparsfix,
                          idparsopt = idparsopt,
                          idparsfix = idparsfix,
                          area = area,
                          obs_richness = obs_richness,
                          num_cycles = num_cycles,
                          f_sar = f_sar)
    if(as.numeric(unlist(out2$fvalues)) > ML) {
      out <- out2
    }
  }
  if(out$conv > 0)
  {
    cat("Optimization has not converged. Try again with different initial values.\n")
    out2 <- NULL
  } else {
    MLtrpars <- as.numeric(unlist(out$par))
    MLpars <- MLtrpars/(1 - MLtrpars)
    MLpars1 <- rep(0,length(idparsopt) + length(idparsfix))
    MLpars1[idparsopt] <- MLpars
    MLpars1[idparsfix] <- parsfix
    ML <- as.numeric(unlist(out$fvalues))
    AIC_c <- AICc(LogLik = ML, k = length(initparsopt), n = length(obs_richness))
    out2 <- c(MLpars1, ML, AIC_c)
    names(out2)[idparsopt] <- names(initparsopt)
    names(out2)[idparsfix] <- names(parsfix)
    names(out2)[(length(out2) - 1):length(out2)] <- c("ML","AICc")
  }
  return(out2)
}

sar_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix, ...)
{
  trpars1 <- rep(0,length(idparsopt) + length(idparsfix))
  trpars1[idparsopt] <- trparsopt
  if(length(idparsfix) != 0)
  {
    trpars1[idparsfix] <- trparsfix
  }
  if(max(trpars1) > 1 || min(trpars1) < 0)
  {
    loglik <- -Inf
  } else {
    pars1 <- trpars1/(1 - trpars1)
    loglik <- sar_loglik(pars = pars1, ...)
    if(is.nan(loglik) || is.na(loglik))
    {
      warning("Parameter values have been used that cause numerical problems.\n")
      loglik <- -Inf
    }
  }
  return(loglik)
}

AIC <- function(LogLik,k){
  aic <- (2 * k) - ( 2 * LogLik)
  return(aic)
}

AICc <- function(LogLik,k,n){
  aicc <- AIC(LogLik,k) + ((2 * k * (k + 1))/(n - k - 1))
  return(aicc)
}

ICweights <- function(IC) {
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights)
}

sar_data <- read.csv('d:/data/ms/lake_area_richness.csv')
area <- sar_data[,4]
obs_richness <- sar_data[,3]
sar_initial_pars <- read.csv('d:/data/ms/sar_initial_pars.csv')
sar_initial_pars21 <- c(c = 8.63, T1 = 4.1, z1 = 1.68, z2 = 5.08)
sar_initial_pars22 <- c(c = 8.63, T1 = 4.1, T2 = 21.14, z1 = 1.68, z2 = 5.08, z3 = 5.08)
sar_initial_pars23 <- c(c = 8.63, T1 = 4.1, z2 = 5.08)
sar_initial_pars24 <- c(c = 8.65, T1 = 4.1, T2 = 21.14, z2 = 5.08, z3 = 5.08)

sar_model_names <- list("f_asymp",
                        "f_betap",
                        "f_chapman",
                        "f_loga",
                        "f_epm1",
                        "f_epm2",
                        "f_gompertz",
                        "f_koba",
                        "f_linear",
                        "f_heleg",
                        "f_monod",
                        "f_mmf",
                        "f_negexpo",
                        "f_p1",
                        "f_p2",
                        "f_power",
                        "f_powerR",
                        "f_ratio",
                        "f_weibull3",
                        "f_weibull4",
                        "f_ContOne",
                        "f_ContTwo",
                        "f_ContOne",
                        "f_ContTwo")
sar_models <- lapply(sar_model_names,get)
out <- list()
for(i in 21:21) {
  if(i <= 20)
  {
    NAs <- which(is.na(sar_initial_pars[i,]))
    if(length(NAs) > 0) {
      initparsopt <- sar_initial_pars[i,-NAs]
    } else {
      initparsopt <- sar_initial_pars[i,]
    }
    idparsopt <- 1:length(initparsopt)
    parsfix <- NULL
    idparsfix <- NULL
  } else {
    if(i == 21) {
      initparsopt <- sar_initial_pars21
      parsfix <- NULL
      idparsopt <- 1:length(initparsopt)
      idparsfix <- NULL
    }
    if(i == 22) {
      initparsopt <- sar_initial_pars22
      parsfix <- NULL
      idparsopt <- 1:length(initparsopt)
      idparsfix <- NULL
    }
    if(i == 23) {
        initparsopt <- sar_initial_pars23
        parsfix <- c(z1 = 0)
        idparsopt <- c(1, 2, 4)
        idparsfix <- 3
    }
    if(i == 24) {
      initparsopt <- sar_initial_pars24
      parsfix <- c(z1 = 0)
      idparsopt <- c(1:3,5:6)
      idparsfix <- 4
    }
  }
  names_pars <- names(initparsopt)
  initparsopt <- as.numeric(initparsopt)
  names(initparsopt) <- names_pars
  names_pars <- names(parsfix)
  parsfix <- as.numeric(parsfix)
  names(parsfix) <- names_pars
  out[[i]] <- list(initpars = rep(0,length(idparsopt) + length(idparsfix)), model = NULL,fit = NULL)
  out[[i]]$initpars[idparsopt] <- initparsopt
  if(length(idparsfix) > 0) out[[i]]$initpars[idparsfix] <- parsfix
  out[[i]]$model <- sar_models[[i]]
  out[[i]]$fit <- sar_ML(f_sar = sar_models[[i]],
                         initparsopt = initparsopt,
                         parsfix = parsfix,
                         idparsopt = idparsopt,
                         idparsfix = idparsfix,
                         area = area,
                         obs_richness = obs_richness,
                         trial_settings = c(n = 10, sd = 0.1))
}
fit_results <- as.data.frame(matrix(NA, nrow = 24, ncol = 10))
names(fit_results) <- c("model", "c","d or T1","f or T2","z or z1","z2","z3","ML","AICc","AICcweights")
for(i in 21:21) {
    vec <- out[[i]]$fit
    lvec <- length(vec)
    fit_results[i,1] <- sar_model_names[[i]]
    fit_results[i,8:9] <- vec[(length(vec) - 1):length(vec)]
    if(!is.na(vec["c"])) fit_results[i,"c"] <- vec["c"]
    if(!is.na(vec["d"])) fit_results[i,"d or T1"] <- vec["d"]
    if(!is.na(vec["f"])) fit_results[i,"f or T2"] <- vec["f"]
    if(!is.na(vec["z"])) fit_results[i,"z or z1"] <- vec["z"]
    if(!is.na(vec["z1"])) fit_results[i,"z or z1"] <- vec["z1"]
    if(!is.na(vec["z2"])) fit_results[i,"z2"] <- vec["z2"]
    if(!is.na(vec["z3"])) fit_results[i,"z3"] <- vec["z3"]
    if(!is.na(vec["T1"])) fit_results[i,"d or T1"] <- vec["T1"]
    if(!is.na(vec["T2"])) fit_results[i,"f or T2"] <- vec["T2"]
}
fit_results[,10] <- ICweights(IC = fit_results[,9])
fit_results[23,1] <- "f_ZslopeOne"
fit_results[24,1] <- "f_ZslopeTwo"
fit_results
new_order <- rev(order(fit_results[,10]))
fit_results[new_order,]
write.csv(fit_results, "d:/data/ms/sar_results.csv", row.names = FALSE)
area2 <- exp(seq(from = log(min(area)), to = log(max(area)), length.out = 500))
pdf("d:/data/ms/sar.pdf", width = 7, height = 5)
par(mfrow = c(1,1))
plot(sort(area2),f_asymp(area = sort(area2),pars = as.numeric(fit_results[1,c(2,3,5)])), xaxt = 'n', type = 'l',log = 'x', col = 'red', lwd = 2, xlab = 'Area',ylab = 'Richness',xlim = c(0.2,600), ylim = c(0,40))
ticks <- axTicks(1)
tick_labels <- c(sprintf("%.1f", ticks[1:2]),sprintf("%.0f", ticks[3:11]))
axis(side = 1, at = ticks, labels = tick_labels)
lines(sort(area2),f_p1(area = sort(area2),pars = as.numeric(fit_results[14,c(2,3,5)])), col = 'blue', lwd = 2)
lines(sort(area2),f_ContOne(area = sort(area2),pars = as.numeric(fit_results[21,c(2,3,5,6)])), col = 'green', lwd = 2)
lines(sort(area2),f_ratio(area = sort(area2),pars = as.numeric(fit_results[18,c(2,3,5)])), col = 'black', lwd = 3)
points(area,obs_richness, pch = 19, cex = 0.5)
legend(x = 0.2,y = 40, legend = c("ratio","ContOne","p1","asymp"), fill = c('black','green','blue','red'), bty = 'n')
dev.off()

