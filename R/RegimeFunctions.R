#' Calibrate the relative bank strength to observed conditions
#' 
#' \code{cal_pattern} is a function that calculates the channel dimensions
#' for a wide range of relative bank strength values, given a set of input
#' parameters (Q, S, D50, D84, D90, fs). The output data can be used to 
#' help calibrate the bank strength parameter used to match the model to an
#' existing river reach
#' 
#' @param Q formative discharge (m3/s)
#' @param S  energy gradient for the stream reach (m/m)
#' @param D50  median size of the bed surface (mm) influencing transport rate
#' @param D84  84th percentile of the bed surface (mmm) influencing resistance
#' @param D90  90th percentile of the bed surface (mm) influencing stability
#' @param Plot a logical value that controls the output of summary plots
#' @param fs  percent of the bed surface that is covered by sand (percent)
#' @export 

cal_pattern = function(Q, S, D50, D84, D90, Plot = T, fs = 0){
  library(lattice)
  mu = seq(1,5,0.1)
  all_runs = est_pattern(Q,mu[1],S,D50,D84,D90,fs)
  
  for (i in seq(2,length(mu))){
    all_runs[i,]= est_pattern(Q,mu[i],S,D50,D84,D90,fs)
  }
  all_runs$mu = mu
  if(Plot == T){
      as.factor(all_runs$N)
      
      print(xyplot(W~mu, data = all_runs, groups = N, 
                   ylab = 'Total wetted width (m)',
                   xlab = 'Relative bank strength (mu)',
                   auto.key = list(space = "top", points = TRUE, columns = 5)
      )
      )
    }
  return(all_runs)
}

#' Run a Monte Carlo simulation predicting stable channel characteristics
#' 
#' \code{sim_pattern} is a function that launches a Monte Carlo simulation to
#' estimate the range of channel conditions that could be expected for a given
#' set of input parameters (Q, mu, S, D50, D84, D90, fs) with a specified level
#' of uncertainty for the input parameters (expressed as a proportion varying 
#' from 0 to 1.0)
#' 
#' @param Q formative discharge (m3/s)
#' @param mu  bank strength, relative to bed (afer Millar, 2005)
#' @param S  energy gradient for the stream reach (m/m)
#' @param D50  median size of the bed surface (mm) influencing transport rate
#' @param D84  84th percentile of the bed surface (mmm) influencing resistance
#' @param D90  90th percentile of the bed surface (mm) influencing stability
#' @param Uncert uncertainty for the inputs, expressed as a proportion between 0 and 1.0
#' @param Nmc number of Monte Carlo simulations to run (default is 100)
#' @param Plot a logical value that controls the output of summary statistics and plots
#' @param fs  percent of the bed surface that is covered by sand (percent)
#' @export 
#' 
sim_pattern = function(Q, mu, S, D50, D84, D90, Uncert = 0.1,  Nmc = 100, Plot = T, fs = 0){
  library(lattice)
  
  Qmc = runif(Nmc, (1-Uncert)*Q, (1+Uncert)*Q)
  Smc = runif(Nmc, (1-Uncert)*S, (1+Uncert)*S)
  D50mc = runif(Nmc, (1-Uncert)*D50, (1+Uncert)*D50)
  D84mc = runif(Nmc, (1-Uncert)*D84, (1+Uncert)*D84)
  D90mc = runif(Nmc, (1-Uncert)*D90, (1+Uncert)*D90)
  
  all_runs = est_pattern(Qmc[1],mu,Smc[1],D50mc[1],D84mc[1],D90mc[1],fs)
  
  for (i in seq(2,Nmc)){
    all_runs[i,]= est_pattern(Qmc[i],mu,Smc[i],D50mc[i],D84mc[i],D90mc[i],fs)
  }
  
  if(Plot == T){
    Nch = sort(unique(all_runs$N))
    
    for (j in Nch){
      print(paste0('NOTE: ', 
                   100*sum(all_runs$N == j)/Nmc, 
                   '% of the simulations had ',
                   j,' anabranches (summary statistics below)'
      )
      )
      
      print(summary(all_runs[all_runs$N == j,]))
      as.factor(all_runs$N)
    }
    print(xyplot(d~W, data = all_runs, groups = N, 
                 xlab = 'Total wetted width (m)',
                 ylab = 'Mean hydraulic depth (m)',
                 auto.key = list(space = "top", points = TRUE, columns = 5)
    )
    )
  }
  return(all_runs)
}
 
#' Run the regime model for stable channel patterns
#' 
#' \code{est_pattern} is a function that runs the UBC Regime Model using
#' the reach average estimates of surface texture, slope, formative discharge
#' and relative bank strength. the function calculates the number of anabranches
#' required for a stable channel pattern, as well as the total channel width (Wt),
#' mean hydraulic depth (d), mean flow velocity (U), and total volumetric transport
#' capacity (Qbt)
#' 
#' @param Q formative discharge (m3/s)
#' @param mu  bank strength, relative to bed (afer Millar, 2005)
#' @param S  energy gradient for the stream reach (m/m)
#' @param D50  median size of the bed surface (mm) influencing transport rate
#' @param D84  84th percentile of the bed surface (mmm) influencing resistance
#' @param D90  90th percentile of the bed surface (mm) influencing stability
#' @param fs  percent of the bed surface that is covered by sand (percent)
#' @export 
#' 
est_pattern = function(Q, mu, S, D50, D84, D90, fs = 0){
  wd_crit = 60
  test = est_regime(Q, mu, S, D50, D84, D90, fs)
  stab_index = (test$W/test$d/wd_crit)
  n = 1
  while(stab_index > 1){
    n = n + 1
    Qan = Q/n
    test = est_regime(Qan, mu, S, D50, D84, D90, fs)
    stab_index = (test$W/test$d/wd_crit)
  }
  solution = test * c(n,1,1,n)
  solution$N = n
  return(solution)
}


#' Run the regime model in basic mode
#' 
#' \code{est_regime} is a function that runs the UBC Regime Model using the
#' reach average estimates of surface texture, slope, formative discharge and
#' relative bank strength. The function calculates the wetted width (W, in m), 
#' mean hydraulic depth (d, in m), mean flow velocity (U, in m/s), and the 
#' sediment transport capacity (Qb, in m3/s)
#' 
#' @param Q formative discharge (m3/s)
#' @param mu  bank strength, relative to bed (afer Millar, 2005)
#' @param S  energy gradient for the stream reach (m/m)
#' @param D50  median size of the bed surface (mm) influencing transport rate
#' @param D84  84th percentile of the bed surface (mmm) influencing resistance
#' @param D90  90th percentile of the bed surface (mm) influencing stability
#' @param fs  percent of the bed surface that is covered by sand (percent)
#' @export 
 
est_regime = function(Q, mu, S, D50, D84, D90, fs = 0){
 
  # guess the convergence tolerance and initial bed width
  Tol = 0.00001
  p = 4 * Q ^ 0.5
  test_plus = find_stable(Q, mu, p*1.001, S, D50, D84, D90, fs)
  test_minus = find_stable(Q, mu, p*0.999, S, D50, D84, D90, fs)
  gradient_1 = test_plus$Qb - test_minus$Qb
  p1 = p
  
  #now move in the direction of the gradient
  if(gradient_1 > 0){p = p + 0.25*p} else {p = p - 0.25*p}
  test_plus = find_stable(Q, mu, p*1.001, S, D50, D84, D90, fs)
  test_minus = find_stable(Q, mu, p*0.999, S, D50, D84, D90, fs)
  gradient_2 = test_plus$Qb - test_minus$Qb
  p2 = p
  
  while(gradient_1/gradient_2 > 0){
    gradient_1 = gradient_2
    p1 = p
    if(gradient_1 > 0){p = p + 0.25*p} else {p = p - 0.25*p}
    test_plus = find_stable(Q, mu, p*1.001, S, D50, D84, D90, fs)
    test_minus = find_stable(Q, mu, p*0.999, S, D50, D84, D90, fs)
    gradient_2 = test_plus$Qb - test_minus$Qb
    p2 = p
  }
  
  p_upper = max(c(p1,p2))
  p_lower = min(c(p1,p2))
  p = 0.5*(p_upper + p_lower)
  converg = (p_upper-p_lower)/p
  
  while(converg > Tol){
    test_plus = find_stable(Q, mu, p*1.001, S, D50, D84, D90, fs)
    test_minus = find_stable(Q, mu, p*0.999, S, D50, D84, D90, fs)
    gradient = test_plus$Qb - test_minus$Qb
    if(gradient > 0){p_lower = p} else {p_upper = p}
    p = 0.5*(p_upper + p_lower)
    converg = (p_upper-p_lower)/p
  }
  solution = find_stable(Q, mu, p, S, D50, D84, D90, fs)
  return(data.frame(solution[c("W", "d", "U","Qb")]))
}

#' Calculate the channel state
#' 
#' \code{channel_state} calculates the characteristics of a trapezoidal channel
#' with specified dimensions and bed material. The calculated characteristics 
#' include wetted width (W), mean depth (d), mean velocity (U), total discharge
#' (Q), shear stress acting on the bed (Tbed), shear stress acting on the bank
#' (Tbank), volumetric transport capacity (Qb)
#' 
#' #' @param Q formative discharge (m3/s)
#' @param p  bed width (m) 
#' @param y maximum flow depth for trapezoid (m)
#' @param b sideslope angle of the banks (degrees)
#' @param S  energy gradient for the stream reach (m/m)
#' @param D50  median size of the bed surface (mm) influencing transport rate
#' @param D84  84th percentile of the bed surface (mmm) influencing resistance
#' @param fs  percent of the bed surface that is covered by sand (percent)
#' 
channel_state = function(p, y, b, S, D50, D84, fs){

# specify the constants and make unit conversions
g <<- 9.81
rho <<- 1000  # water density
Gs <<- 1.65   # submerged specific gravity
a1 = 6.5
a2 = 2.5
b = b * pi / 180
D50 = D50 / 1000
D84 = D84 / 1000
fs = fs / 100 

# calculate the area of flow
dW = y / tan(b)
dP = y / sin(b)
A = y * (p + dW)
R = A / (p + 2 * dP)  # area / perimeter

# use Ferguson 2007 to calculate the stream velocity
Res = a1 * a2 * (R/D84) / (a1^2 + a2^2 * (R/D84)^(5/3))^(1/2)
# use the Keulegan Equation
# Res = (1/0.4)*log(12.2*R/(D84))
U = Res * (g * R * S)^(1/2)

# use the equations from Knight and others to partition stress
SFbank = 10^(-1.4026 * log10(p/(2 * dP) + 1.5) + 0.247)
bed_stress =  g * rho * y * S * (1 - SFbank) * ((p + 2 * dW) / (2 * p) + 0.5)        
bank_stress =  g * rho * y * S * SFbank * (2 * p + 2 * dW)*sin(b)/(4*y);  

# use Wilcock and Crowe to estimate the sediment transport rate
tau_star_ref = 0.021 + 0.015 * exp (-20 * fs)
tau_ref = tau_star_ref * g * Gs * rho * D50
X = bed_stress / tau_ref
if(X < 1.35){
  W_star = 0.002 * X^(7.5)
}else{
  W_star = 14 * (1 - (0.894/(X^0.5)))^(4.5)
}
Qb = p * (W_star / (1.65 * g)) * (bed_stress/rho)^(3/2)
  
#create a list to store data
state = list()
state$W = p + 2 * dW
state$d = A / state$W   # area / width
state$U = U
state$Q = A * U
state$Tbed = bed_stress
state$Tbank = bank_stress
state$Qb = Qb
state$b = b * 180 /pi  # save the value of sideslope and bed width
state$p = p
return(state)
}

#' Solve for water depth
#' 
#' \code{find_Q} uses a midpoint solution approach to find water depth that
#' produces the user-specified stream discharge in a trapezoidal channel with
#' a known geometry.
#' @param Q discharge contained by the channel (m3/s)
#' @param p  bed width (m) 
#' @param b sideslope angle of the banks (degrees)
#' @param S  energy gradient for the stream reach (m/m) 
#' @param D50  median size of the bed surface (mm) influencing transport rate
#' @param D84 84th percentile of the bed surface grain size distribution (mm)
#' @param fs  percent of the bed surface that is covered by sand (percent)
#' 
find_Q = function(Q, p, b, S, D50, D84, fs){
  # find the appropriate depth to satisfy continuity with the
  # specified design flow, Q
  
  # guess the initial channel depth, set up variables
  y = 0.3 * Q ^ 0.3
  deltaX = 0.001 * Q ^ 0.3
  tol = 0.00001
  
  # first test
  test = channel_state(p, y, b, S, D50, D84, fs)
  converg = (test$Q - Q) / Q
  
  while(abs(converg) > tol){
    test2 = channel_state(p,y+deltaX, b, S, D50, D84, fs)
    M = (test2$Q - test$Q)/deltaX
    B = (test$Q - Q) - M*y
    y = -B/M
    test = channel_state(p, y, b, S, D50, D84, fs)
    converg = (test$Q - Q) / Q 
  }
  return(y)
}

#' Solve for the stable channel geometry
#' 
#' \code{find_stable} solves for the stable bank angle for a given channel 
#' width, slope, sediment texture, relative bank strength and discharge. The
#' function returns all of the descriptor variables in \code{channel_state}
#' 
#' @param Q formative discharge (m3/s)
#' @param mu  bank strength, relative to bed (afer Millar, 2005)
#' @param p  bed width (m) 
#' @param S  energy gradient for the stream reach (m/m)
#' @param D50  median size of the bed surface (mm) influencing transport rate
#' @param D84  84th percentile of the bed surface (mmm) influencing resistance
#' @param D90  90th percentile of the bed surface (mm) influencing stability
#' @param fs  percent of the bed surface that is covered by sand (percent)

find_stable = function(Q, mu, p, S, D50, D84, D90, fs){
  # find the stable channel shape for the specified Q and 
  # relative bank strength, mu
 
  # specify constants and set mu to an equivalent mod_phi
  phi <<- 40
  mod_phi <<- atan(mu * tan(phi * pi / 180)) * 180/pi
  Tol = 0.00001
  deltaX = 0.00001*mod_phi
  #tau_star <<- 0.02
  tau_star <<- 0.02 #default value for use with the D90
  
  # set the upper and lower angle limits
  b_upper = mod_phi - deltaX
  b_lower = deltaX
  b = 0.67 * mod_phi
  
  # calculate the bank stability index
  y = find_Q(Q,p,b,S,D50,D84,fs)
  test = channel_state(p,y,b,S,D50,D84,fs)
  bank_crit = g*rho*Gs*(D90/1000) * tau_star * 
              (tan(mod_phi*pi/180) / tan(phi*pi/180)) *
              (1 - (sin(b*pi/180)^2 / sin(mod_phi*pi/180)^2))^0.5
  converg = (test$Tbank - bank_crit) / bank_crit
  
  while(abs(converg) > Tol){
    if(converg > 0){b_upper = b} else {b_lower = b}
    b = 0.5*(b_upper + b_lower)
    y = find_Q(Q,p,b,S,D50,D84,fs)
    test = channel_state(p,y,b,S,D50,D84,fs)
    bank_crit = g*rho*Gs*(D90/1000) * tau_star * 
                (tan(mod_phi*pi/180) / tan(phi*pi/180)) *
                (1 - (sin(b*pi/180)^2 / sin(mod_phi*pi/180)^2))^0.5
    converg = (test$Tbank - bank_crit) / bank_crit
  }
  return(test)
}
