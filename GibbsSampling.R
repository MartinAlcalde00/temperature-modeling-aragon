### IMPLEMENTACIÓN DE GIBBS SAMPLING. MODELO NO JERÁRQUICO.

library(extraDistr)
library(xtable)
library(coda)

# Datos de temperaturas y altitudes

temp <- readRDS("meanTempAragonJJA19562015.rds.ds")
elev <- readRDS("elev.rds")

# Matriz de datos. Filas por años y columnas por localidades. 
# Generamos la fila 1 de y, que se corresponde con y_0,s.
# En general, la fila t refiere a las temperaturas del año t-1.

y <- matrix(nrow = 61, ncol = 18)
y[2:61,] <- temp
y[1,] <- apply(temp, MARGIN = 2, FUN = mean)

# Definimos la matriz que contiene los 60 años (1,...,60) en cada columna.

t <- matrix(data = 1:60, nrow = 60, ncol = 18, byrow = FALSE)

# Definimos la matriz alt que contiene las 18 de las localidades en cada fila.

alt <- matrix(data = elev$altitude, nrow = 60, ncol = 18, byrow = TRUE)

# Parámetros

numberOfSamples <- 200000

beta_0 <- rep(0, numberOfSamples)
alpha <- rep(0, numberOfSamples)
gamma <- rep(0, numberOfSamples)
rho <- rep(0, numberOfSamples)
sigma2 <- rep(0, numberOfSamples)

beta_02 <- rep(0, numberOfSamples)
alpha2 <- rep(0, numberOfSamples)
gamma2 <- rep(0, numberOfSamples)
rho2 <- rep(0, numberOfSamples)
sigma22 <- rep(0, numberOfSamples)

# Valores iniciales de los parámetros.

beta_0[1] <- 0
alpha[1] <- 0
gamma[1] <- 0
sigma2[1] <- 1
rho[1] <- 0

beta_02[1] <- 1
alpha2[1] <- 1
gamma2[1] <- 1
sigma22[1] <- 2
rho2[1] <- 0.2

# Iteraciones de Gibbs sampling.

sumSquares <- sum(y[ -61,]^2) # Esta constante es necesaria la función sampleRho
sumHeightsSq <- sum(alt[1,]^2)*60 # Esta constante es necesaria la función sampleGamma

for (i in 2:numberOfSamples) {
  beta_0[i] <- sampleBeta0(0, 5000, alpha[i-1], gamma[i-1], rho[i-1], 
                             sigma2[i-1]) 
  alpha[i] <- sampleAlpha(0, 5000, beta_0[i], gamma[i-1], rho[i-1], 
                            sigma2[i-1])
  gamma[i] <- sampleGamma(0, 5000, beta_0[i], alpha[i], rho[i-1], 
                            sigma2[i-1])
  rho[i] <- sampleRho(beta_0[i], alpha[i], gamma[i], sigma2[i-1])
    
  sigma2[i] <- sampleSigma2(1, 1, beta_0[i], alpha[i], rho[i], 
                              gamma[i])
}

for (i in 2:numberOfSamples) {
  beta_02[i] <- sampleBeta0(0, 5000, alpha2[i-1], gamma2[i-1], rho2[i-1],
                            sigma22[i-1])
  alpha2[i] <- sampleAlpha(0, 5000, beta_02[i], gamma2[i-1], rho2[i-1], 
                           sigma22[i-1])
  gamma2[i] <- sampleGamma(0, 5000, beta_02[i], alpha2[i], rho2[i-1], 
                           sigma22[i-1])
  rho2[i] <- sampleRho(beta_02[i], alpha2[i], gamma2[i], sigma22[i-1])
  
  sigma22[i] <- sampleSigma2(1, 1, beta_02[i], alpha2[i], rho2[i], 
                             gamma2[i])
}




# Implementación de las distribuciones condicionales a posteriori y muestreo.

sampleBeta0 <- function(mu, tau2, alpha, gamma, rho, sigma2) {
  suma <- sum(y[-1,] - alpha * t - rho * y[-61,] - gamma * alt) / sigma2
  precision <- (1 / tau2 + 1080 / sigma2)
  media <- (mu / tau2 + suma) / precision
  return(rnorm(n = 1, mean = media, sd = 1 / sqrt(precision)))
}

sampleAlpha <- function(mu, tau2, beta_0, gamma, rho, sigma2) {
  suma <- sum(t * (y[-1,] - beta_0 - rho * y[-61,] - gamma * alt)) / sigma2
  precision <- (1 / tau2 +  1328580 / sigma2)
  media <- (mu / tau2 + suma) / precision
  return(rnorm(n = 1, mean = media, sd = 1 / sqrt(precision)))
}

sampleGamma <- function(mu, tau2, beta_0, alpha, rho, sigma2) {
  suma <- sum(alt * (y[-1,] - beta_0 - alpha * t - rho * y[-61,])) / sigma2
  precision <- (1 / tau2 + sumHeightsSq / sigma2)
  media <- (mu / tau2 + suma) / precision
  return(rnorm(n = 1, mean = media, sd = 1 / sqrt(precision)))
}

sampleRho <- function(beta_0, alpha, gamma, sigma2) {
  suma <- sum(y[-61,] * (y[-1,] - beta_0 - alpha * t - gamma * alt))
  media <- suma / sumSquares
  varianza <- sigma2 / sumSquares
  return(rtnorm(n = 1, mean = media, sd = sqrt(varianza), -1, 1))
}

sampleSigma2 <- function(n, s2, beta_0, alpha, rho, gamma) {
  suma <- sum((y[-1,] - beta_0 - alpha * t - rho * y[-61,] - gamma * alt)^2)
  scale2 <- n * s2 + suma
  return(scale2 / rchisq(n = 1, df = 1080 + n))
}

# Traceplots de ambas cadenas

plot(beta_0[1:1], type='l', col = 'blue', xlab = 'Iteraciones', 
     ylab = expression(beta[0]))
lines(beta_02[100001:200000], col = 'red')
plot(alpha[100001:200000], type='l', col = 'blue', xlab = 'Iteraciones', 
     ylab = expression(alpha))
lines(alpha2[100001:200000], col = 'red')
plot(rho[100001:200000], type='l', col = 'blue', xlab = 'Iteraciones', 
     ylab = expression(rho))
lines(rho2[100001:200000], col = 'red')
plot(gamma[100001:200000], type='l', col = 'blue', xlab = 'Iteraciones', 
     ylab = expression(gamma))
lines(gamma2[100001:200000], col = 'red')
plot(sigma2[100001:200000], type='l', col = 'blue', xlab = 'Iteraciones', 
     ylab = expression(sigma^2))
lines(sigma22[100001:200000], col = 'red')

# Factor de reducción de escala pontecial RHat.

gelman.diag(mcmc.list(as.mcmc(beta_0[100001:200000]), 
                      as.mcmc(beta_02[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
gelman.diag(mcmc.list(as.mcmc(alpha[100001:200000]), 
                      as.mcmc(alpha2[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
gelman.diag(mcmc.list(as.mcmc(rho[100001:200000]), 
                      as.mcmc(rho2[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
gelman.diag(mcmc.list(as.mcmc(gamma[100001:200000]), 
                      as.mcmc(gamma2[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
gelman.diag(mcmc.list(as.mcmc(sigma2[100001:200000]), 
                      as.mcmc(sigma22[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)


# Densidades marginales a posteriori

plot(density(beta_0[100001:200000]), col = 'green')
lines(density(beta_02[100001:200000]), col = 'blue')
plot(density(alpha[100001:200000]), col = 'green')
lines(density(alpha2[100001:200000]), col = 'blue')
plot(density(rho[100001:200000]), col = 'green')
lines(density(rho2[100001:200000]), col = 'blue')
plot(density(gamma[100001:200000]), col = 'green')
lines(density(gamma2[100001:200000]), col = 'blue')
plot(density(sigma2[100001:200000]), col = 'green')
lines(density(sigma22[100001:200000]), col = 'blue')

# Valores aproximados de la esperanza a posteriori

mean(beta_0[100001:200000])
mean(alpha[100001:200000])
mean(rho[100001:200000]) 
mean(gamma[100001:200000])
mean(sigma2[100001:200000])

#Intervalos de credibilidad

matrix <- cbind(beta_0[100001:200000], alpha[100001:200000], 
                rho[100001:200000], gamma[100001:200000],
                sigma2[100001:200000])

lowerBoundsCI <- apply(matrix, MARGIN = 2, FUN = function(x) 
  quantile(x, probs = 0.025)) 

upperBoundsCI <- apply(matrix, MARGIN = 2, FUN = function(x) 
  quantile(x, probs = 0.975))

#Resumen numérico completo

posteriorMeans <-c(mean(beta_0[100001:200000]),
                   mean(alpha[100001:200000]),
                   mean(rho[100001:200000]), 
                   mean(gamma[100001:200000]),
                   mean(sigma2[100001:200000]))

summaryTable <- cbind(posteriorMeans, lowerBoundsCI, upperBoundsCI)

xtable(summaryTable, caption = 'Resumen de los resultados', digits = c(4,4,4,4))

 #Guardando las gráficas
pdf(file = 'beta_0Densidades.pdf')
plot(density(beta_0[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~beta[0]))
lines(density(beta_02[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'alphaDensidades.pdf')
plot(density(alpha[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~alpha))
lines(density(alpha2[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'rhoDensidades.pdf')
plot(density(rho[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~rho))
lines(density(rho2[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'gammaTDensidades.pdf')
plot(density(gamma[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~gamma))
lines(density(gamma2[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'sigma2TDensidades.pdf')
plot(density(sigma2[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~sigma^2))
lines(density(sigma22[100001:200000]), col = 'blue')
dev.off()

###

pdf(file = 'rhoTP.pdf')
plot(rho[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(rho))
lines(rho2[100001:200000], col = 'blue')
dev.off()

pdf(file = 'gammaTP.pdf')
plot(gamma[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(gamma))
lines(gamma2[100001:200000], col = 'blue')
dev.off()

pdf(file = 'sigma2TP.pdf')
plot(sigma2[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(sigma^2))
lines(sigma22[100001:200000], col = 'blue')
dev.off()

pdf(file = 'beta_0TP.pdf')
plot(beta_0[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(beta[0]))
lines(beta_02[100001:200000], col = 'blue')
dev.off()

pdf(file = 'alphaTP.pdf')
plot(alpha[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(alpha))
lines(alpha2[100001:200000], col = 'blue')
dev.off()

cosas <- cbind(gelman.diag(mcmc.list(as.mcmc(beta_0[100001:200000]), 
                                     as.mcmc(beta_02[100001:200000])), 
                           confidence = 0.95, autoburnin = FALSE),
               gelman.diag(mcmc.list(as.mcmc(alpha[100001:200000]), 
                                     as.mcmc(alpha2[100001:200000])), 
                           confidence = 0.95, autoburnin = FALSE),
               gelman.diag(mcmc.list(as.mcmc(rho[100001:200000]), 
                                     as.mcmc(rho2[100001:200000])), 
                           confidence = 0.95, autoburnin = FALSE),
               gelman.diag(mcmc.list(as.mcmc(gamma[100001:200000]), 
                                     as.mcmc(gamma2[100001:200000])), 
                           confidence = 0.95, autoburnin = FALSE),
               gelman.diag(mcmc.list(as.mcmc(sigma2[100001:200000]), 
                                     as.mcmc(sigma22[100001:200000])), 
                           confidence = 0.95, autoburnin = FALSE))