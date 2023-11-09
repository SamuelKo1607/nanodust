remove(list=ls())
#setwd(dir = "C:\\Users\\skoci\\Documents\\nanodust")
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(INLA)
library(hexbin)
require(hexbin)
require(lattice)
#require(RColorBrewer)
#library(vioplot)




###################################
### DEFINE MODEL with R-generic ###
###################################


three_component_model = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                 "log.prior", "quit"), theta=NULLL){

  envir <-parent.env(environment())
  prec.high = exp(15)
  
  
  rate <- function(v_sc_r, v_sc_t, r_sc, v_sc_x, v_sc_y, v_sc_z,  #covariates
                   l_bg, l_isd, l_b, v_b_r, e_v, e_r,      #hyperparameters
                   S=8,                                    #sc parameters
                   v_isd=26, phi_isd=259, th_isd=8,        #isd parameters
                   v_b_a=9, v_earth_a=29.8){               #beta met. parameters
    
    deg2rad <- function(deg){
      rad = deg/180*pi
      return(rad)
    }
    
    #background
    L_bg = l_bg
    
    #interstellar contribution
    v_isd_x = -v_isd * sin(deg2rad(90-th_isd)) * cos(deg2rad(phi_isd))
    v_isd_y = -v_isd * sin(deg2rad(90-th_isd)) * sin(deg2rad(phi_isd))
    v_isd_z = -v_isd * cos(deg2rad(90-th_isd)) 
    L_isd = l_isd * (   ( ( v_isd_x - v_sc_x )^2 
                          +  ( v_isd_y - v_sc_y )^2 
                          +  ( v_isd_z - v_sc_z )^2 )^0.5 
    ) / (v_isd)
    
    #beta meteoroid contribution
    ksi = -2 - e_r
    r_factor = r_sc/1
    v_factor = ( (
      ( v_sc_r - ( v_b_r*(r_factor^ksi)  ) )^2 
      + ( v_sc_t - ( v_b_a*(r_factor^(-1)) ) )^2
    )^0.5 
    ) / ( (
      ( v_b_r )^2 
      + ( v_earth_a - v_b_a )^2
    )^0.5 
    )
    L_b = l_b * (v_factor)^e_v * (r_factor)^e_r 
    
    #normalization to hourly rate, while L_i are in m^-2 s^-1
    hourly_rate = 3600 * S * ( L_bg + L_isd + L_b )
    return(hourly_rate)
  }

  
  interpret.theta <- function(){
    return(list(l_bg  = exp(theta[1L]), 
                l_isd = exp(theta[2L]),
                l_b   = exp(theta[3L]),
                v_b_r = theta[4L],
                e_v   = theta[5L], 
                e_r   = theta[6L]
               ))
  }
  
  graph <-function(){
    G <- Diagonal(n = length(vt), x=1)
    return(G)
  }
  
  Q <- function(){
    #prec.high <- interpret.theta()$prec
    Q <- prec.high*graph()
    return(Q)
  }
  
  mu <- function(){
    par = interpret.theta()
    return(log( rate(#covariates
                     vr, vt, r, vx, vy, vz, 
                     #hyperparameters
                     par$l_bg, par$l_isd, par$l_b, par$v_b_r, par$ e_v, par$e_r     
                     )
              ))
  }
  
  log.norm.const <-function(){
    return(numeric(0))
  }
  
  # Log-prior for thetas
  log.prior <- function(){
    par = interpret.theta()
    
    #nice priors
    val <- (dgamma(par$l_bg,  shape = 2,    scale = 1e-6, log=TRUE) + theta[1L]+
            dgamma(par$l_isd, shape = 2,    scale = 1e-5, log=TRUE) + theta[2L]+
            dgamma(par$l_b,   shape = 2,    scale = 1e-4, log=TRUE) + theta[3L]+
            dnorm(par$v_b_r,  mean  = 60,    sd   = 5,     log=TRUE) +
            dnorm(par$e_v,    mean  = 2.2,   sd   = 0.05,  log=TRUE) +
            dnorm(par$e_r,    mean  = -1.65, sd   = 0.05,  log=TRUE) 
           )
  
    return(val)
  }
  
  # Initial values of theta
  initial <- function(){
    return(c(log(1e-6),
             log(1e-5),
             log(1e-4),
             60,
             2.2,
             -1.65
            )
          )   
  }
  
  quit <-function(){
    return(invisible())
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return(val)
}






###################################
########## Load the data ##########
###################################

mydata = read.csv(file = 'data_synced\\flux_readable.csv')
names(mydata)[c(2,3,4,5,6,9,10,11)] = c("flux",
                                        "vr",
                                        "vt",
                                        "r",
                                        "exposure",
                                        "vx",
                                        "vy",
                                        "vz")
n = length(mydata$vr)
mydata$idx = 1:n 









###################################
########## Run the model ##########
###################################


rgen = inla.rgeneric.define(model = three_component_model, 
                            vr = mydata$vr, 
                            vt = mydata$vt, 
                            r  = mydata$r,
                            vx = mydata$vx,
                            vy = mydata$vy,
                            vz = mydata$vz)
result = inla(flux ~ -1 + f(idx, model = rgen),
              data = mydata, family = "poisson", E = exposure, 
              control.compute = list(cpo=TRUE, dic=TRUE, config = TRUE),
              safe = TRUE, verbose = TRUE)

summary(result)

hist(result$cpo$pit)     # ok 
result$cpo$failure       # also OK
pit = result$cpo$pit
save(pit, file = "998_generated\\inla\\pit.RData")
  
plot(mydata$flux/mydata$exposure, ylab="counts/E")
lines(result$summary.fitted.values$mean, col=2, lwd=3)
  
# Posterior means of the hyperparameters
inla.emarginal(function(x) exp(x), result$marginals.hyperpar$`Theta1 for idx`)
inla.emarginal(function(x) exp(x), result$marginals.hyperpar$`Theta2 for idx`)
inla.emarginal(function(x) exp(x), result$marginals.hyperpar$`Theta3 for idx`)
result$summary.hyperpar$mean[4]
result$summary.hyperpar$mean[5]
result$summary.hyperpar$mean[6]


plot(exp(result$marginals.hyperpar$`Theta1 for idx`[1:43]),result$marginals.hyperpar$`Theta1 for idx`[44:86],main = "l_bg")
plot(exp(result$marginals.hyperpar$`Theta2 for idx`[1:43]),result$marginals.hyperpar$`Theta2 for idx`[44:86],main = "l_isd")
plot(exp(result$marginals.hyperpar$`Theta3 for idx`[1:43]),result$marginals.hyperpar$`Theta3 for idx`[44:86],main = "l_b")
plot((result$marginals.hyperpar$`Theta4 for idx`[1:43]),result$marginals.hyperpar$`Theta4 for idx`[44:86],main = "v_b_r")
plot((result$marginals.hyperpar$`Theta5 for idx`[1:43]),result$marginals.hyperpar$`Theta5 for idx`[44:86],main = "e_v")
plot((result$marginals.hyperpar$`Theta6 for idx`[1:43]),result$marginals.hyperpar$`Theta6 for idx`[44:86],main = "e_r")






###################################
###### Sample the posterior  ######
###################################


s = inla.hyperpar.sample(1000000, result)

l_bg  = exp(s[,1])
l_isd = exp(s[,2])
l_b   = exp(s[,3])
v_b_r =     s[,4]
e_v   =     s[,5]
e_r   =     s[,6]


hexbinplot(v_b_r~l_isd, 
           data=data.frame(v_b_r,l_isd), 
           colramp=colorRampPalette(c("grey", "yellow")),
           main="joint histogram v_b_r, l_isd" ,  
           xlab="l_isd", 
           ylab="v_b_r" ,
           panel=function(x, y, ...)
           {
             panel.hexbinplot(x, y, ...)
             panel.abline(v=c(mean(l_isd)), h=c(mean(v_b_r)), col="black", lwd=2, lty=3)
           }
)

save(l_bg, 
     l_isd, 
     l_b, 
     v_b_r, 
     e_v, 
     e_r, 
     file = "998_generated\\inla\\sample.RData")






