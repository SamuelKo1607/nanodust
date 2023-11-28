remove(list=ls())
setwd(dir = "C:\\Users\\skoci\\Documents\\nanodust")
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# install.packages("distr")
# install.packages('dplyr')

library(INLA)
library(hexbin)
require(hexbin)
require(lattice)
library(distr)
library(dplyr)
#require(RColorBrewer)
#library(vioplot)



###################################
### DEFINE MODEL with R-generic ###
###################################

three_component_model <- function(cmd = c("graph", "Q", "mu", "initial", 
                                          "log.norm.const", "log.prior", "quit",
                                          "rate", 
                                          "prior.l_bg",
                                          "prior.l_isd",
                                          "prior.l_b",
                                          "prior.v_b_r",
                                          "prior.e_v",
                                          "prior.e_r"), 
                                  theta=NULL, feed_x=NULL){

  envir <- parent.env(environment())
  prec.high = exp(15)
  
  prior.l_bg <- function(l_bg=feed_x){
    return(dgamma(l_bg,  shape = 2,    scale = 1e-5, log=TRUE))
  }
  prior.l_isd <- function(l_isd=feed_x){
    return(dgamma(l_isd, shape = 2,    scale = 1e-5, log=TRUE))
  }
  prior.l_b <- function(l_b=feed_x){
    return(dgamma(l_b,   shape = 2,    scale = 1e-4, log=TRUE))
  }
  prior.v_b_r <- function(v_b_r=feed_x){
    return(dnorm(v_b_r,  mean  = 60,    sd   = 5,    log=TRUE))
  }
  prior.e_v <- function(e_v=feed_x){
    return(dnorm(e_v,    mean  = 2.2,   sd   = 0.05, log=TRUE))
  }
  prior.e_r <- function(e_r=feed_x){
    return(dnorm(e_r,    mean  = -1.65, sd   = 0.05, log=TRUE))
  }
  
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
    val <- (prior.l_bg(  par$l_bg)    + theta[1L] +
            prior.l_isd( par$l_isd)   + theta[2L] +
            prior.l_b(   par$l_b)     + theta[3L] +
            prior.v_b_r( par$v_b_r)   +
            prior.e_v(   par$e_v)     + 
            prior.e_r(   par$e_r)
           )
  
    return(val)
  }
  
  # Initial values of theta
  initial <- function(){
    #initial values set to the maxima a priori
    return(c(log(optimize(prior.l_bg, interval = c(0, 1e-2), maximum = TRUE, tol=1e-9)$maximum),
             log(optimize(prior.l_isd, interval = c(0, 1e-2), maximum = TRUE, tol=1e-9)$maximum),
             log(optimize(prior.l_b, interval = c(0, 1e-2), maximum = TRUE, tol=1e-9)$maximum),
             optimize(prior.v_b_r, interval = c(0, 1000), maximum = TRUE, tol=1e-6)$maximum,
             optimize(prior.e_v, interval = c(-100, 100), maximum = TRUE, tol=1e-6)$maximum,
             optimize(prior.e_r, interval = c(-100, 100), maximum = TRUE, tol=1e-6)$maximum
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

# for inner only fit
# mydata = filter(mydata, r <= 0.8)

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
#save(pit, file = "998_generated\\inla\\pit.RData")

#plotting
par(mfrow = c(1, 1))
plot(mydata$flux/mydata$exposure, ylab="counts/E")
lines(result$summary.fitted.values$mean, col=2, lwd=3)
lines(30+mydata$flux/mydata$exposure-result$summary.fitted.values$mean, col="blue")

span = round(max(abs(mydata$flux/mydata$exposure-result$summary.fitted.values$mean),40)+0.5)
hist(mydata$flux/mydata$exposure-result$summary.fitted.values$mean,
     breaks=c(-span:span),
     main="")
mtext(paste("residuals histogram, stdev = ",
            as.character(sqrt(var(mydata$flux/mydata$exposure-result$summary.fitted.values$mean))),
            ", log(mlik) = ",
            result$mlik[1]), side=3)
  
# Posterior means of the hyperparameters
inla.emarginal(function(x) exp(x), result$marginals.hyperpar$`Theta1 for idx`)
inla.emarginal(function(x) exp(x), result$marginals.hyperpar$`Theta2 for idx`)
inla.emarginal(function(x) exp(x), result$marginals.hyperpar$`Theta3 for idx`)
result$summary.hyperpar$mean[4]
result$summary.hyperpar$mean[5]
result$summary.hyperpar$mean[6]

# Create a layout with one column and six rows
par(mfrow = c(6, 1), mar = c(2, 2, 1, 1))
# Plot each function in a separate row
plot(exp(result$marginals.hyperpar$`Theta1 for idx`[1:43]),result$marginals.hyperpar$`Theta1 for idx`[44:86])
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]), y = par("usr")[4] - 0.3 * diff(par("usr")[3:4]), labels = "l_bg", col = "red", cex = 1.5)
plot(exp(result$marginals.hyperpar$`Theta2 for idx`[1:43]),result$marginals.hyperpar$`Theta2 for idx`[44:86])
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]), y = par("usr")[4] - 0.3 * diff(par("usr")[3:4]), labels = "l_isd", col = "red", cex = 1.5)
plot(exp(result$marginals.hyperpar$`Theta3 for idx`[1:43]),result$marginals.hyperpar$`Theta3 for idx`[44:86])
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]), y = par("usr")[4] - 0.3 * diff(par("usr")[3:4]), labels = "l_b", col = "red", cex = 1.5)
plot((result$marginals.hyperpar$`Theta4 for idx`[1:43]),result$marginals.hyperpar$`Theta4 for idx`[44:86])
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]), y = par("usr")[4] - 0.3 * diff(par("usr")[3:4]), labels = "v_b_r", col = "red", cex = 1.5)
plot((result$marginals.hyperpar$`Theta5 for idx`[1:43]),result$marginals.hyperpar$`Theta5 for idx`[44:86])
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]), y = par("usr")[4] - 0.3 * diff(par("usr")[3:4]), labels = "e_v", col = "red", cex = 1.5)
plot((result$marginals.hyperpar$`Theta6 for idx`[1:43]),result$marginals.hyperpar$`Theta6 for idx`[44:86])
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]), y = par("usr")[4] - 0.3 * diff(par("usr")[3:4]), labels = "e_r", col = "red", cex = 1.5)
# Reset the layout to the default (1x1)
par(mfrow = c(1, 1))



###################################
### Priors and post. evaluation ###
###################################

#posteriors to be saved in X/Y form
#l_bg
fx_l_bg = exp(result$marginals.hyperpar$`Theta1 for idx`[1:43])
fy_l_bg = result$marginals.hyperpar$`Theta1 for idx`[44:86]
#l_isd
fx_l_isd = exp(result$marginals.hyperpar$`Theta2 for idx`[1:43])
fy_l_isd = result$marginals.hyperpar$`Theta2 for idx`[44:86]
#l_b
fx_l_b = exp(result$marginals.hyperpar$`Theta3 for idx`[1:43])
fy_l_b = result$marginals.hyperpar$`Theta3 for idx`[44:86]
#v_b_r
fx_v_b_r = result$marginals.hyperpar$`Theta4 for idx`[1:43]
fy_v_b_r = result$marginals.hyperpar$`Theta4 for idx`[44:86]
#e_v
fx_e_v = result$marginals.hyperpar$`Theta5 for idx`[1:43]
fy_e_v = result$marginals.hyperpar$`Theta5 for idx`[44:86]
#e_r
fx_e_r = result$marginals.hyperpar$`Theta6 for idx`[1:43]
fy_e_r = result$marginals.hyperpar$`Theta6 for idx`[44:86]

#priors to be saved in X/Y form
#log priors extracted
prior.l_bg  <- function(x){ 
  return((three_component_model(cmd="prior.l_bg", feed_x=x))) }
prior.l_isd <- function(x){
  return((three_component_model(cmd="prior.l_isd",feed_x=x))) }
prior.l_b   <- function(x){
  return((three_component_model(cmd="prior.l_b",  feed_x=x))) }
prior.v_b_r <- function(x){
  return((three_component_model(cmd="prior.v_b_r",feed_x=x))) }
prior.e_v   <- function(x){
  return((three_component_model(cmd="prior.e_v",  feed_x=x))) }
prior.e_r   <- function(x){
  return((three_component_model(cmd="prior.e_r",  feed_x=x))) }
#priors evaluated
#l_bg
p_l_bg_max = optimize(prior.l_bg, interval = c(0, 1), maximum = TRUE, tol=1e-9)$maximum
px_l_bg = seq(p_l_bg_max/1000, p_l_bg_max*100, length.out = 100000)
py_l_bg = exp(prior.l_bg(px_l_bg))
#l_isd
p_l_isd_max = optimize(prior.l_isd, interval = c(0, 1), maximum = TRUE, tol=1e-9)$maximum
px_l_isd = seq(p_l_isd_max/1000, p_l_isd_max*100, length.out = 100000)
py_l_isd = exp(prior.l_isd(px_l_isd))
#l_b
p_l_b_max = optimize(prior.l_b, interval = c(0, 1), maximum = TRUE, tol=1e-9)$maximum
px_l_b = seq(p_l_b_max/1000, p_l_b_max*100, length.out = 100000)
py_l_b = exp(prior.l_b(px_l_b))
#v_b_r
p_l_v_b_r_max = optimize(prior.v_b_r, interval = c(0, 10000), maximum = TRUE, tol=1e-9)$maximum
px_v_b_r = seq(p_l_v_b_r_max/100, p_l_v_b_r_max*10, length.out = 10000)
py_v_b_r = exp(prior.v_b_r(px_v_b_r))
#e_v
p_e_v_max = optimize(prior.e_v, interval = c(0, 10), maximum = TRUE, tol=1e-9)$maximum
px_e_v = seq(p_e_v_max/5, p_e_v_max*5, length.out = 10000)
py_e_v = exp(prior.e_v(px_e_v))
#e_r
p_e_r_max = optimize(prior.e_r, interval = c(-10, 0), maximum = TRUE, tol=1e-9)$maximum
px_e_r = seq(p_e_r_max*5, p_e_r_max/5, length.out = 10000)
py_e_r = exp(prior.e_r(px_e_r))



###################################
######## Sample posterior  ########
###################################

s = inla.hyperpar.sample(1000000, result)

sample_l_bg  = exp(s[,1])
sample_l_isd = exp(s[,2])
sample_l_b   = exp(s[,3])
sample_v_b_r =     s[,4]
sample_e_v   =     s[,5]
sample_e_r   =     s[,6]

hexbinplot(sample_v_b_r~sample_l_isd, 
           data=data.frame(sample_v_b_r,sample_l_isd), 
           colramp=colorRampPalette(c("grey", "yellow")),
           main="joint histogram v_b_r, l_isd" ,  
           xlab="l_isd", 
           ylab="v_b_r" ,
           panel=function(x, y, ...)
           {
             panel.hexbinplot(x, y, ...)
             panel.abline(v=c(mean(sample_l_isd)), h=c(mean(sample_v_b_r)), col="black", lwd=2, lty=3)
           }
)

hexbinplot(sample_l_bg~sample_l_isd, 
           data=data.frame(sample_l_bg,sample_l_isd), 
           colramp=colorRampPalette(c("grey", "blue")),
           main="joint histogram l_bg, l_isd" ,  
           xlab="l_isd", 
           ylab="l_bg" ,
           panel=function(x, y, ...)
           {
             panel.hexbinplot(x, y, ...)
             panel.abline(v=c(mean(sample_l_isd)), h=c(mean(sample_l_bg)), col="black", lwd=2, lty=3)
           }
)



###################################
#### Save everything relevant #####
###################################

model_definition <- deparse(three_component_model)

current_time <- Sys.time()
formatted_time <- format(current_time, "%Y%m%d%H%M%S")

save(sample_l_bg,       #sampled posterior
     sample_l_isd, 
     sample_l_b, 
     sample_v_b_r, 
     sample_e_v, 
     sample_e_r,
     fx_l_bg,           #evaluated posterior
     fy_l_bg,
     fx_l_isd,
     fy_l_isd,
     fx_l_b,
     fy_l_b,
     fx_v_b_r,
     fy_v_b_r,
     fx_e_v,
     fy_e_v,
     fx_e_r,
     fy_e_r,
     px_l_bg,           #evaluated prior
     py_l_bg,
     px_l_isd,
     py_l_isd,
     px_l_b,
     py_l_b,
     px_v_b_r,
     py_v_b_r,
     px_e_v,
     py_e_v,
     px_e_r,
     py_e_r,
     model_definition,  #the model definition
     file = paste("998_generated\\inla\\sample_",formatted_time,".RData",
                  sep = ""))






