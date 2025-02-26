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
