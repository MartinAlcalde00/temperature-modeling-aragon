### IMPLEMENTACIÓN DE GIBBS SAMPLING. MODELO NO JERÁRQUICO.

library(extraDistr)
library(xtable)
library(coda)
library(rockchalk)
# Semilla

set.seed(780646)

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

tau2 <- rep(0, numberOfSamples)
beta_0 <- rep(0, numberOfSamples)
gamma <- rep(0, numberOfSamples)
rho <- rep(0, numberOfSamples)
sigma2 <- rep(0, numberOfSamples)
# La primera fila de deltas corresponde a delta_0. En general, la fila t 
# corresponde a delta_{t-1}
deltas <- matrix(data = 0, nrow = numberOfSamples, ncol = 60)

tau22 <- rep(0, numberOfSamples)
beta_02 <- rep(0, numberOfSamples)
gamma2 <- rep(0, numberOfSamples)
rho2 <- rep(0, numberOfSamples)
sigma22 <- rep(0, numberOfSamples)
# La primera fila de deltas2 corresponde a delta_0. En general, la fila t 
# corresponde a delta_{t-1}
deltas2 <- matrix(data = 0, nrow = numberOfSamples, ncol = 60)

# Valores iniciales de los parámetros.

tau2[1] <- 1
beta_0[1] <- 0
gamma[1] <- 0
sigma2[1] <- 1
rho[1] <- 0
deltas[1,] <- rep(0, 60)

tau22[1] <- 2
beta_02[1] <- 1
gamma2[1] <- 1
sigma22[1] <- 2
rho2[1] <- 0.2
deltas2[1,] <- rep(1, 60)

# Iteraciones de Gibbs sampling.

sumSquares <- sum(y[ -61,]^2) # Esta constante es necesaria la función sampleRho
sumHeightsSq <- sum(alt[1,]^2)*60 # Esta constante es necesaria la función sampleGamma

for (i in 2:numberOfSamples) {
  for (j in 1:60) {
    deltas[i,j] <- sampleDeltas(j, beta_0[i-1], tau2[i-1], gamma[i-1], 
                                rho[i-1], sigma2[i-1])
    mDeltas <- matrix(data = deltas[i,], nrow = 60, ncol = 18)
  }                          
  gamma[i] <- sampleGamma(0, 5000, mDeltas, rho[i-1], sigma2[i-1])
  
  rho[i] <- sampleRho(mDeltas, gamma[i], sigma2[i-1])
  
  sigma2[i] <- sampleSigma2(1, 1, mDeltas, rho[i], gamma[i])
  
  beta_0[i] <- sampleBeta_0(deltas[i,], tau2[i-1])
  
  tau2[i] <- sampleTau2(1, 1, beta_0[i], deltas[i,])
}

for (i in 2:numberOfSamples) {
  for (j in 1:60) {
    deltas2[i,j] <- sampleDeltas(j, beta_02[i-1], tau22[i-1], gamma2[i-1],
                                 rho2[i-1], sigma22[i-1])
    mDeltas <- matrix(data = deltas2[i,], nrow = 60, ncol = 18)
  } 
  gamma2[i] <- sampleGamma(0, 5000, mDeltas, rho2[i-1], 
                          sigma22[i-1])
  rho2[i] <- sampleRho(mDeltas, gamma2[i], sigma22[i-1])
  
  sigma22[i] <- sampleSigma2(1, 1, mDeltas, rho2[i], gamma2[i])
  
  beta_02[i] <- sampleBeta_0(deltas2[i,], tau22[i-1])
  
  tau22[i] <- sampleTau2(1, 1, beta_02[i], deltas2[i,])
}

# Implementación de las distribuciones condicionales a posteriori y muestreo.

sampleDeltas <- function(j, beta_0, tau2, gamma, rho, sigma2) {
  suma <- sum((y[j + 1,] - rho * y[j,] - gamma * alt[1,]) / sigma2)
  precision <- (1 / tau2 + 18 / sigma2)
  media <- (beta_0 / tau2 + suma) / precision
  return(rnorm(n = 1, mean = media, sd = 1 / sqrt(precision)))
}

sampleGamma <- function(mu, tau2, deltas, rho, sigma2) {
  suma <- sum(alt * (y[-1,] - deltas - rho * y[-61,])) / sigma2
  precision <- (1 / tau2 + sumHeightsSq / sigma2)
  media <- (mu / tau2 + suma) / precision
  return(rnorm(n = 1, mean = media, sd = 1 / sqrt(precision)))
}

sampleRho <- function(deltas, gamma, sigma2) {
  suma <- sum(y[-61,] * (y[-1,] - deltas - gamma * alt))
  media <- suma / sumSquares
  varianza <- sigma2 / sumSquares
  return(rtnorm(n = 1, mean = media, sd = sqrt(varianza), -1, 1))
}

sampleSigma2 <- function(n, s2, deltas, rho, gamma) {
  suma <- sum((y[-1,] - deltas - rho * y[-61,] - gamma * alt)^2)
  scale2 <- n * s2 + suma
  return(scale2 / rchisq(n = 1, df = 1080 + n))
}

sampleBeta_0 <- function(deltas, tau2) {
  return(rnorm(n = 1, mean = (sum(deltas)) / 60, sd = sqrt(tau2 / 60)))
}

sampleTau2 <- function(n, s2, beta_0, deltas) {
  scale <- sum((deltas - beta_0)^2) + s2
  return(scale / rchisq(n = 1, df = n + 60))
}

# Traceplots de ambas cadenas

plot(rho[100001:200000], type='l', col = 'blue', xlab = 'Iteraciones', 
     ylab = expression(rho))
lines(rho2[100001:200000], col = 'green')
plot(gamma[100001:200000], type='l', col = 'green', xlab = 'Iteraciones', 
     ylab = expression(gamma))
lines(gamma2[100001:200000], col = 'blue')
plot(sigma2[100001:200000], type='l', col = 'green', xlab = 'Iteraciones', 
     ylab = expression(sigma^2))
lines(sigma22[100001:200000], col = 'blue')
plot(beta_0[100001:200000], type='l', col = 'green', xlab = 'Iteraciones', 
     ylab = expression(beta[0]))
lines(beta_02[100001:200000], col = 'blue')
plot(tau2[100001:200000], type='l', col = 'green', xlab = 'Iteraciones', 
     ylab = expression(tau^2))
lines(tau22[100001:200000], col = 'blue')

# Factor de reducción de escala pontecial RHat.

gelman.diag(mcmc.list(as.mcmc(rho[100001:200000]), 
                      as.mcmc(rho2[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
gelman.diag(mcmc.list(as.mcmc(gamma[100001:200000]), 
                      as.mcmc(gamma2[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
gelman.diag(mcmc.list(as.mcmc(sigma2[100001:200000]), 
                      as.mcmc(sigma22[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
gelman.diag(mcmc.list(as.mcmc(beta_0[100001:200000]), 
                      as.mcmc(beta_02[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
gelman.diag(mcmc.list(as.mcmc(tau2[100001:200000]), 
                      as.mcmc(tau22[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE)
for (i in 1:60) {
  gelman.diag(mcmc.list(as.mcmc(deltas[100001:200000, i]), 
                        as.mcmc(deltas2[100001:200000, i])), 
              confidence = 0.95, autoburnin = FALSE)
}

# Densidades marginales a posteriori

plot(density(rho[100001:200000]), col = 'blue')
lines(density(rho2[100001:200000]), col = 'red')
plot(density(gamma[100001:200000]), col = 'blue')
lines(density(gamma2[100001:200000]), col = 'red')
plot(density(sigma2[100001:200000]), col = 'blue')
lines(density(sigma22[100001:200000]), col = 'red')
plot(density(beta_0[100001:200000]), col = 'blue')
lines(density(beta_02[100001:200000]), col = 'red')
plot(density(tau2[100001:200000]), col = 'blue')
lines(density(tau22[100001:200000]), col = 'red')

# Boxplot de los deltas

boxplot(deltas[100001:200000, ], outline = FALSE)

# Estimaciones de las esperanzas a posteriori

mean(rho[100001:200000]) 
mean(gamma[100001:200000])
mean(sigma2[100001:200000])
mean(beta_0[100001:200000])
mean(tau2[100001:200000])

# Intervalos de credibilidad

matrix <- cbind(rho[100001:200000], gamma[100001:200000],
                sigma2[100001:200000], beta_0[100001:200000],
                tau2[100001:200000])

lowerBoundsCI <- apply(matrix, MARGIN = 2, FUN = function(x) 
  quantile(x, probs = 0.025)) 

upperBoundsCI <- apply(matrix, MARGIN = 2, FUN = function(x) 
  quantile(x, probs = 0.975))

### MÁS COSAS

posteriorMeans <- c(mean(rho[100001:200000]), 
                    mean(gamma[100001:200000]),
                    mean(sigma2[100001:200000]),
                    mean(beta_0[100001:200000]),
                    mean(tau2[100001:200000]))

summaryTable <- cbind(posteriorMeans, lowerBoundsCI, upperBoundsCI)

xtable(summaryTable, caption = 'Resumen de los resultados', digits = c(4,4,4,4))

#Guardando las gráficas

pdf(file = 'rhoDensidades2.pdf')
plot(density(rho[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~rho))
lines(density(rho2[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'gammaTDensidades2.pdf')
plot(density(gamma[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~gamma))
lines(density(gamma2[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'sigma2TDensidades2.pdf')
plot(density(sigma2[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~sigma^2))
lines(density(sigma22[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'beta_0Densidades2.pdf')
plot(density(beta_0[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~beta[0]))
lines(density(beta_02[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'tau2Densidades2.pdf')
plot(density(tau2[100001:200000]), col = 'green', xlab = NA, ylab = NA, 
     main = bquote('Densidad'~tau^2))
lines(density(tau22[100001:200000]), col = 'blue')
dev.off()

pdf(file = 'deltas.pdf')
boxplot(deltas[100001:200000, ], outline = FALSE)
dev.off()



cosas <- cbind(gelman.diag(mcmc.list(as.mcmc(rho[100001:200000]), 
                      as.mcmc(rho2[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE),
gelman.diag(mcmc.list(as.mcmc(gamma[100001:200000]), 
                      as.mcmc(gamma2[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE),
gelman.diag(mcmc.list(as.mcmc(sigma2[100001:200000]), 
                      as.mcmc(sigma22[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE),
gelman.diag(mcmc.list(as.mcmc(beta_0[100001:200000]), 
                      as.mcmc(beta_02[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE),
gelman.diag(mcmc.list(as.mcmc(tau2[100001:200000]), 
                      as.mcmc(tau22[100001:200000])), 
            confidence = 0.95, autoburnin = FALSE))

cosas2 <- NULL

for (i in 1:60) {
  cosas2 <- cbind(cosas2, )
}

pdf(file = 'rhoTP2.pdf')
plot(rho[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(rho))
lines(rho2[100001:200000], col = 'blue')
dev.off()

pdf(file = 'gammaTP2.pdf')
plot(gamma[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(gamma))
lines(gamma2[100001:200000], col = 'blue')
dev.off()

pdf(file = 'sigma2TP2.pdf')
plot(sigma2[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(sigma^2))
lines(sigma22[100001:200000], col = 'blue')
dev.off()

pdf(file = 'beta_0TP2.pdf')
plot(beta_0[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(beta[0]))
lines(beta_02[100001:200000], col = 'blue')
dev.off()

pdf(file = 'tau2TP2.pdf')
plot(tau2[100001:200000], type = 'l', col = 'green', 
     xlab = 'Iteraciones', ylab = expression(tau^2))
lines(tau22[100001:200000], col = 'blue')
dev.off()


