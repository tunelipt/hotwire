


hwSignalConditioner <- function(offset=0, gain=1, lowPassFilter=30000){

  config <- list(offser=offset, gain=gain, lowPassFilter=lowPassFilter)
  class(config) <- "SignalConditioner"
  return(config)
}


hwBridgeConfig <- function(shape="Film", bridgeRatio="1:20", ampGain=8, filter=3, cableCompensation=0, coolingInterval=2){

  config <- list(shape=shape, bridgeRatio=bridgeRatio, ampGain=ampGain, filter=filter, cableCompensation=0, coolingInterval=coolingInterval)
  class(config) <- "BridgeConfig"

  return(config)
}



hwCable <- function(name="A1863", L=4, R = 0.2, impedence=50,  connectors=2, description=""){

  cable <- list(name=name, L=L, R=R, impedence=impedence, connectors=connectors, description=description)
  class(cable) <- "Cable"
  return(cable)
}


hwSupport <- function(name='55H21', D=4, L=235, Lc=765, contacts=1, R=0.44){
  support <- list(name=name, D=D, L=L, Lc=Lc, contacts=contacts, R=R)
  class(support) <- "Support"
}

hwWire <- function(R=3.4, alpha=0.4, T0=25, RL=0.5, overheat=0.8, X=pi/2, Y=0, Z=0, k=0.04, h=1, Tmax=300, Tamax=150){


  wire <- list(R=R, alpha=alpha, T0=T0, RL=RL, overheat=overheat, X=X, Y=Y, Z=Z, k=k, h=h, Tmax=Tmax, Tamax=Tamax)
  class(wire) <- "Wire"
  return(wire)
}


hwProbe <- function(support, cable, bridge, signal, ...){

  sensors <- list(...)

  n <- length(sensors)


  probe <- list(sensors=sensors, nsensors=n, support=support, cable=cable, bridge=bridge, signal=signal)
  if (n==1)
    class(probe) <- "Probe1D"
  else if (n==2)
    class(probe) <- "Probe2D"
  else if (n==3)
    class(probe) <- "Probe3D"
  else
    stop("Sensor deve ter entre 1 e 3 fios!")
  

  return(probe)
}

hwRefTemp <- function(probe)
  return(probe$sensors[[1]]$T0)

hwOperatingTemp <- function(probe, sensor=1){
  wire <- probe$sensors[[sensor]]

  R0 <- wire$R
  T0 <- wire$T0
  alpha <- wire$alpha
  a <- wire$overheat
  Rw <- R0*(a+1)
  Tw <- (Rw-R0) / (alpha/100*R0) + T0
  return(Tw)
}


  

hwCalibrProbe1D <- function(probe, E, Ucal, temperature=NULL, curvefit=kingVel){
  if (is.null(temperature)){
    temperature <- hwRefTemp(probe)
  }
  T0 <- hwRefTemp(probe)
  Tw <- hwOperatingTemp(probe, 1)
  
  calibr <- hwCalibr(E, Ucal, temperature, T0, Tw,  curvefit)
  return(calibr)
}

hwCalibrProbe3D <- function(probe, E1, Ucal1, E2, Ucal2, E3, Ucal3, temperature=NULL, curvefit=kingVel){
  if (is.null(temperature)){
    temperature <- hwRefTemp(probe)
  }
  T0 <- hwRefTemp(probe)
  Tw1 <- hwOperatingTemp(probe, 1)
  Tw2 <- hwOperatingTemp(probe, 2) 
  Tw3 <- hwOperatingTemp(probe, 3)
  calibr1 <- hwCalibr(E1, Ucal1, temperature, T0, Tw1,  curvefit)
  calibr2 <- hwCalibr(E2, Ucal2, temperature, T0, Tw2,  curvefit)
  calibr3 <- hwCalibr(E3, Ucal3, temperature, T0, Tw3,  curvefit)
  k1 <- probe$sensors[[1]]$k
  k2 <- probe$sensors[[2]]$k
  k3 <- probe$sensors[[3]]$k
  h1 <- probe$sensors[[1]]$h
  h2 <- probe$sensors[[2]]$h
  h3 <- probe$sensors[[3]]$h

  a1 <- sapply(probe$sensors[[1]][c('X', 'Y', 'Z')], cos)
  a2 <- sapply(probe$sensors[[2]][c('X', 'Y', 'Z')], cos)
  a3 <- sapply(probe$sensors[[3]][c('X', 'Y', 'Z')], cos)
  
  f <- 1/( k1*(k2*k3-h3) - h2*k3 + h1*(h2*h3-k2) + 1 )
  fe1 <- (k1*a1[1]^2 +    a2[1]^2 + h1*a3[1]^2)
  fe2 <- (h2*a1[1]^2 + k2*a2[1]^2 +    a3[1]^2)
  fe3 <- (   a1[1]^2 + h3*a2[1]^2 + k3*a3[1]^2)
  
  function(E1, E2, E3, Ta=T0){
    Ue1 <- calibr1(E1, Ta)^2 * fe1
    Ue2 <- calibr2(E2, Ta)^2 * fe2
    Ue3 <- calibr3(E3, Ta)^2 * fe3
    
    U1 <- sqrt(f * abs( (k2*k3-h3)* Ue1 + (h1*h3-k3)*Ue2 + (1-h1*k2)*Ue3 ) ) 
    U2 <- sqrt(f * abs( (1-h2*k3)* Ue1 + (k1*k3-h1)*Ue2 + (h1*h2 - k1)*Ue3 ) )
    U3 <- sqrt(f * abs( (h2*h3-k2)* Ue1 + (1-h3*k1)*Ue2 + (k1*k2-h2)*Ue3 ) )


    return(cbind(-a1[1]*U1 - a2[1]*U2 - a3[1]*U3,
                 -a1[2]*U1 - a2[2]*U2 - a3[2]*U3,
                 -a1[3]*U1 - a2[3]*U2 - a3[3]*U3))
  }
}

hwCalibrProbe2D <- function(probe, E1, Ucal1, E2, Ucal2, temperature=NULL, curvefit=kingVel){
  if (is.null(temperature)){
    temperature <- hwRefTemp(probe)
  }
  T0 <- hwRefTemp(probe)
  Tw1 <- hwOperatingTemp(probe, 1)
  Tw2 <- hwOperatingTemp(probe, 2)
  calibr1 <- hwCalibr(E1, Ucal1, temperature, T0, Tw1,  curvefit)
  calibr2 <- hwCalibr(E2, Ucal2, temperature, T0, Tw2,  curvefit)
  k1 <- probe$sensors[[1]]$k
  k2 <- probe$sensors[[2]]$k
  a1 <- probe$sensors[[1]]$X
  a2 <- probe$sensors[[2]]$X
  f <- 1/(k1*k2-1)
  fe1 <- cos(a1)^2 * (1+k1)
  fe2 <- cos(a2)^2 * (1+k2)
  function(E1, E2, Ta=T0){
    Ue1 <- calibr1(E1, Ta)^2*fe1
    Ue2 <- calibr2(E2, Ta)^2*fe2
    U1 <- sqrt(f * (k2*Ue1 - Ue2))
    U2 <- sqrt(f * (-Ue1 +k1*Ue2))

    return(cbind(U1*cos(a1) - U2*cos(a2), U1*sin(a1) - U2*sin(a2)))
  }


}


 

