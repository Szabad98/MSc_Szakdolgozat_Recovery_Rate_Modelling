source("/MCMC_fuggvenyek.R")

x1 <- rnorm(24,0,1)
x2 <- rnorm(24,0,1)

library(forecast)
library(tseries)

acf(x1)
pacf(x1)
auto.arima(x1,max.d = 0, max.order = 10)
arma_1 <- arima(x1,order= c(1,0,0),include.mean=FALSE)
beta_1 <- 1
checkresiduals(arma_1)
u1 <- as.numeric(residuals(arma_1))

acf(x2)
pacf(x2)
auto.arima(x2,max.d = 0, max.order = 10)
arma_2 <- arima(x2,order= c(1,0,0),include.mean=FALSE)
beta_2 <- 1
checkresiduals(arma_2)
u2 <- as.numeric(residuals(arma_2))


PD_0 <- 91202/25572087
y_D <- -qnorm(PD_0,0,1)

#kell a három függvény a nagy fájlból

eredmeny <- MCMC_fv(u1,u2,x1,x2,beta_1,beta_2,y_D,
                    omega_start = 0.15,
                    rho1_start = 0.15,
                    rho2_start = 0.15,
                    epszilon1 = 0.002,
                    epszilon2 = 0.002,
                    epszilon3 = 0.002 ,
                    futas = 5000,
                    pozitivitas = FALSE,
                    meredekseg = 1.362159)

hossz <- nrow(eredmeny)

omega_abrazol(eredmeny,hossz)
rho_1_abrazol(eredmeny,hossz)
rho_2_abrazol(eredmeny,hossz)


# ilyen nagy epszilon értékek mellett nem működik a modell
# kisebb epszilonokra is 0-ba száll rho1 és rho2
# ezek alapján Witzany kezdeti paraméterezése sem működik
# ezálal Witzany modellje nem mindig helytálló
