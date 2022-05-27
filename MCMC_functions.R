# functions for MCMC sim #


#likelihood function
likelihood <- function(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1,rho2,meredekseg){ 
  # az likelihood mindig egyről induljon
  lh = 1
  for (i in 1:24) {
    szamlalo1 = dnorm(u1[i])
    nevezo1 = beta1 * sqrt(abs(rho1)/(1-abs(rho1))) * dnorm((sqrt(abs(rho1)) * x1[i] - y_D) / sqrt(1-rho1))
    szamlalo2 =  dnorm(u2[i])
    nevezo2 = sqrt(1-omega^2) * beta2 * sqrt(abs(rho2)) * meredekseg * (1/sqrt(2*pi)) *
              exp((-1*rho2*((omega*x1[i]+sqrt(1-omega^2)*x2[i])^2))/(2*(2-rho2)))
    
    lh = lh * (szamlalo1/nevezo1) * (szamlalo2/nevezo2)
  }
  
  lh
} 


# MCMC sim function
MCMC_fv <- function(u1,u2,x1,x2,beta1,beta2,y_D,omega_start,rho1_start,rho2_start,
                    epszilon1,epszilon2,epszilon3,epszilon4,futas,pozitivitas,meredekseg){
  
  if (pozitivitas==TRUE) {
    hatar1 <- 0
    hatar2 <- 0
    random_beta <- rbeta(futas,1,1)
    random_beta2 <- 2*rbeta(futas,1,1)-1
  } else {
    hatar1 <- -1
    hatar2 <- -2
    random_beta <- 2*rbeta(futas,1,1)-1
    random_beta2 <- 2*rbeta(futas,1,1)-1
  }
  
  # Create an empty data.frame to store the accepted parameter values for each iteration.  
  # Remember: "the posterior probability is just an updated version of the prior"  
  posterior <- data.frame()  
  
  # Set the starting value of p  
  omega <- omega_start
  rho1 <- rho1_start
  rho2 <- rho2_start
  
  # Start the loop (MCMC)  
  for (i in 1:futas) {  
    # Obtain a new proposal value for omega  
    omega_prime <- omega + runif(1, -epszilon1,epszilon1)
    # Avoid values out of the range -1 - 1  
    if (omega_prime < -1) {
      omega_prime <- -2 - omega_prime
    }  
    if (omega_prime > 1) {
      omega_prime <- 2 - omega_prime
    }
    # Compute the acceptance probability using our likelihood function and the  
    # beta(1,1) distribution as our prior probability.  
    R1 <- likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega_prime,rho1,rho2,meredekseg) /
      likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1,rho2,meredekseg) * 
      dbeta((omega_prime+1)/2,1,1)/dbeta((omega+1)/2,1,1) 
    # Accept or reject the new value of omega  
    if (R1 > 1) {
      R1 <- 1
    }  
    random1 <- runif (1,0,1)  
    if (random1 < R1) {  
      omega <- omega_prime  
    }
    # Store the likelihood of the accepted parameter and its value  
    posterior[i,1] <- log(likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1,rho2,meredekseg))  
    posterior[i,2] <- omega
    
    
    # STEP2
    # Obtain a new proposal value for rho1  
    rho1_prime <- rho1 + runif(1, -epszilon2,epszilon2)
    # Avoid values out of the range -1 - 1 
    if (rho1_prime < hatar1) {
      rho1_prime <- hatar2 - rho1_prime
    }  
    if (rho1_prime >= 1) {
      rho1_prime <- 1.9999 - rho1_prime
    }
    # Compute the acceptance probability using our likelihood function and the  
    # beta(1,1) distribution as our prior probability. 
    R2 <- likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1_prime,rho2,meredekseg) /
      likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1,rho2,meredekseg) * 
      dbeta((rho1_prime+1)/2,1,1)/dbeta((rho1+1)/2,1,1)
    # Accept or reject the new value of rho1  
    if (R2 > 1) {
      R2 <- 1
    }  
    random2 <- runif (1,0,1)  
    if (random2 < R2) {  
      rho1 <- rho1_prime 
    }
    # Store the likelihood of the accepted parameter and its value  
    posterior[i,3] <- log(likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1,rho2,meredekseg))  
    posterior[i,4] <- rho1
    
    
    # STEP3
    # Obtain a new proposal value for rho2  
    rho2_prime <- rho2 + runif(1, -epszilon3,epszilon3)
    # Avoid values out of the range -1 - 1 
    if (rho2_prime < hatar1) {
      rho2_prime <- hatar2 - rho2_prime
    }  
    if (rho2_prime >= 1) {
      rho2_prime <- 1.9999 - rho2_prime
    }
    # Compute the acceptance probability using our likelihood function and the  
    # beta(1,1) distribution as our prior probability. 
    R3 <- likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1,rho2_prime,meredekseg) /
      likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1,rho2,meredekseg) * 
      dbeta((rho2_prime+1)/2,1,1)/dbeta((rho2+1)/2,1,1)
    # Accept or reject the new value of rho2  
    if (R3 > 1) {
      R3 <- 1
    }  
    random3 <- runif (1,0,1)  
    if (random3 < R3) {  
      rho2 <- rho2_prime 
    }
    # Store the likelihood of the accepted parameter and its value  
    posterior[i,5] <- log(likelihood(u1,u2,x1,x2,beta1,beta2,y_D,omega,rho1,rho2,meredekseg))  
    posterior[i,6] <- rho2
    
    posterior[i,7] <- random_beta[i]
    posterior[i,8] <- random_beta2[i]
  } 
  
  # get the posterior df back  
  posterior  
}


# MCMC parameterbecsles abrazolasa (3 db fuggveny)
library(ggplot2)
library(patchwork)

omega_abrazol <- function(eredmeny, hossz){
  #omega
  plot1 <- ggplot(eredmeny, aes(1:hossz, V2)) +
    # geom_point()+
    geom_line(lwd=1.5) +
    geom_line(y=mean(eredmeny$V2), col="red",lwd=1.5)+
    labs(
      title = "omega ertekenek valtozasa",
      x = "szimulacio",
      y = "omega"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14L,face = "bold",hjust = 0.5),
      axis.title.y = element_text(size = 12L,hjust = 0.5),
      axis.title.x = element_text(size = 12L,hjust = 0.5)
    ) 
  
  plot2 <- ggplot(eredmeny) +
    geom_density(aes(V2),lwd=1.5, col="red",adjust=2.5) +
    geom_density(aes(V8),lwd=1.5, col="blue",linetype="dotted",adjust=2) +
    labs(
      title = "prior vs posterior surusegfv.",
      x = "omega",
      y = "suruseg"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14L,face = "bold",hjust = 0.5),
      axis.title.y = element_text(size = 12L,hjust = 0.5),
      axis.title.x = element_text(size = 12L,hjust = 0.5)
    ) +
    xlim(min(min(eredmeny$V2),min((eredmeny$V7))), max(max(eredmeny$V2),max((eredmeny$V7))))
  # xlim(-1,1)
  
  plot1+plot2
}


rho_1_abrazol <- function(eredmeny, hossz){
  #rho1
  plot3 <- ggplot(eredmeny, aes(1:hossz, V4)) +
    # geom_point()+
    geom_line(lwd=1.5) +
    geom_line(y=mean(eredmeny$V4), col="red",lwd=1.5)+
    labs(
      title = "rho1 ertekenek valtozasa",
      x = "szimulacio",
      y = "rho1"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14L,face = "bold",hjust = 0.5),
      axis.title.y = element_text(size = 12L,hjust = 0.5),
      axis.title.x = element_text(size = 12L,hjust = 0.5)
    ) 
  
  plot4 <- ggplot(eredmeny) +
    geom_density(aes(V4),lwd=1.5, col="red",adjust=2.5) +
    geom_density(aes(V7),lwd=1.5, col="blue",linetype="dotted",adjust=2) +
    labs(
      title = "prior vs posterior surusegfv.",
      x = "rho1",
      y = "suruseg"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14L,face = "bold",hjust = 0.5),
      axis.title.y = element_text(size = 12L,hjust = 0.5),
      axis.title.x = element_text(size = 12L,hjust = 0.5)
    ) +
    xlim(min(min(eredmeny$V4),min((eredmeny$V7))), max(max(eredmeny$V4),max((eredmeny$V7))))
  # xlim(-1,1)
  
  plot3+plot4
}

rho_2_abrazol <- function(eredmeny, hossz){
  #rho2
  plot5 <- ggplot(eredmeny, aes(1:hossz, V6)) +
    # geom_point()+
    geom_line(lwd=1.5) +
    geom_line(y=mean(eredmeny$V6), col="red",lwd=1.5)+
    labs(
      title = "rho2 ertekenek valtozasa",
      x = "szimulacio",
      y = "rho2"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14L,face = "bold",hjust = 0.5),
      axis.title.y = element_text(size = 12L,hjust = 0.5),
      axis.title.x = element_text(size = 12L,hjust = 0.5)
    ) 
  
  plot6 <- ggplot(eredmeny) +
    geom_density(aes(V6),lwd=1.5, col="red",adjust=2.5) +
    geom_density(aes(V7),lwd=1.5, col="blue",linetype="dotted",adjust=2) +
    labs(
      title = "prior vs posterior surusegfv.",
      x = "rho2",
      y = "suruseg"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14L,face = "bold",hjust = 0.5),
      axis.title.y = element_text(size = 12L,hjust = 0.5),
      axis.title.x = element_text(size = 12L,hjust = 0.5)
    ) +
    xlim(min(min(eredmeny$V6),min((eredmeny$V7))), max(max(eredmeny$V6),max((eredmeny$V7))))
  # xlim(-1,1)
  
  plot5+plot6
}
