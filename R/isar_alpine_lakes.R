isar_alpine_lakes <- function()
{
  isar_data <- read.csv('d:/data/ms/lake_area_richness.csv')
  area <- isar_data[,4]
  obs_richness <- isar_data[,3]
  isar_initial_pars <- read.csv('d:/data/ms/isar_initial_pars.csv')
  isar_initial_pars21 <- c(c = 8.63, T1 = 4.1, z1 = 1.68, z2 = 5.08)
  isar_initial_pars22 <- c(c = 8.63, T1 = 4.1, T2 = 21.14, z1 = 1.68, z2 = 5.08, z3 = 5.08)
  isar_initial_pars23 <- c(c = 8.63, T1 = 4.1, z2 = 5.08)
  isar_initial_pars24 <- c(c = 8.65, T1 = 4.1, T2 = 21.14, z2 = 5.08, z3 = 5.08)

  isar_model_names <- list("f_asymp",
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
  isar_models <- lapply(isar_model_names,get)
  out <- list()
  for(i in 21:21) {
    if(i <= 20)
    {
      NAs <- which(is.na(isar_initial_pars[i,]))
      if(length(NAs) > 0) {
        initparsopt <- isar_initial_pars[i,-NAs]
      } else {
        initparsopt <- isar_initial_pars[i,]
      }
      idparsopt <- 1:length(initparsopt)
      parsfix <- NULL
      idparsfix <- NULL
    } else {
      if(i == 21) {
        initparsopt <- isar_initial_pars21
        parsfix <- NULL
        idparsopt <- 1:length(initparsopt)
        idparsfix <- NULL
      }
      if(i == 22) {
        initparsopt <- isar_initial_pars22
        parsfix <- NULL
        idparsopt <- 1:length(initparsopt)
        idparsfix <- NULL
      }
      if(i == 23) {
        initparsopt <- isar_initial_pars23
        parsfix <- c(z1 = 0)
        idparsopt <- c(1, 2, 4)
        idparsfix <- 3
      }
      if(i == 24) {
        initparsopt <- isar_initial_pars24
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
    out[[i]]$model <- isar_models[[i]]
    out[[i]]$fit <- isar_ML(f_isar = isar_models[[i]],
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
    fit_results[i,1] <- isar_model_names[[i]]
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
  write.csv(fit_results, "d:/data/ms/isar_results.csv", row.names = FALSE)
  area2 <- exp(seq(from = log(min(area)), to = log(max(area)), length.out = 500))
  pdf("d:/data/ms/isar.pdf", width = 7, height = 5)
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
}
