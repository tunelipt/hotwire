

#' Configuração do condicionador de sinais
#' @export
hwSignalConditioner <- function(offset=0, gain=1, lowPassFilter=30000){

  config <- list(offser=offset, gain=gain, lowPassFilter=lowPassFilter)
  class(config) <- "SignalConditioner"
  return(config)
}


#' Configuração da ponte
#' @export
hwBridgeConfig <- function(shape="Film", bridgeRatio="1:20", ampGain=8, filter=3, cableCompensation=0, coolingInterval=2){

  config <- list(shape=shape, bridgeRatio=bridgeRatio, ampGain=ampGain, filter=filter, cableCompensation=0, coolingInterval=coolingInterval)
  class(config) <- "BridgeConfig"

  return(config)
}



#' Configuração do cabo
#' @export
hwCable <- function(name="A1863", L=4, R = 0.2, impedence=50,  connectors=2, description=""){

  cable <- list(name=name, L=L, R=R, impedence=impedence, connectors=connectors, description=description)
  class(cable) <- "Cable"
  return(cable)
}


#' Configuração do suporte
#' @export
hwSupport <- function(name='55H21', D=4, L=235, Lc=765, contacts=1, R=0.44){
  support <- list(name=name, D=D, L=L, Lc=Lc, contacts=contacts, R=R)
  class(support) <- "Support"
}

#' Parâmetros de configuração de um fio.
#' @export
hwWire <- function(R=3.4, alpha=0.4, T0=25, RL=0.5, overheat=0.8, X=pi/2, Y=0, Z=0, k=0.04, h=1, Tmax=300, Tamax=150){


  wire <- list(R=R, alpha=alpha, T0=T0, RL=RL, overheat=overheat, X=X, Y=Y, Z=Z, k=k, h=h, Tmax=Tmax, Tamax=Tamax)
  class(wire) <- "Wire"
  return(wire)
}

#' Configuração de uma sonda.
#' @export
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


#' Calcula a temperatura de referência.
#' @export
hwRefTemp <- function(probe)
  return(probe$sensors[[1]]$T0)

#' Calcula a temperatura de operação.
#' @export
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


  
#' Calibração de sensor 1D
#' @export
hwCalibrProbe1D <- function(probe, E, Ucal, temperature=NULL, curvefit=kingVel){
  if (is.null(temperature)){
    temperature <- hwRefTemp(probe)
  }
  T0 <- hwRefTemp(probe)
  Tw <- hwOperatingTemp(probe, 1)
  
  calibr <- hwCalibr(E, Ucal, temperature, T0, Tw,  curvefit)
  return(calibr)
}

#' Calibração de sensor 3D
#' @export
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
  A <- matrix(c(k1, 1, h1, h2, k2, 1, 1, h3, k3), 3, 3, byrow=TRUE)
  A <- solve(A)
  aa <- cos(d2r(54.74))^2
  fe1 <- (k1 +    1 +  h1)*aa
  fe2 <- (h2 +    k2 + 1)*aa
  fe3 <- (1  +    h3 + k3)*aa
  
  #fe1 <- k1*a1[1]^2 +    a2[1]^2 + h1*a3[1]^2
  #fe2 <- h2*a1[1]^2 + k2*a2[1]^2 +    a3[1]^2
  #fe3 <-    a1[1]^2 + h3*a2[1]^2 + k3*a3[1]^2
  function(E1, E2, E3, Ta=T0){
    Ue1 <- calibr1(E1, Ta)^2 * fe1
    Ue2 <- calibr2(E2, Ta)^2 * fe2
    Ue3 <- calibr3(E3, Ta)^2 * fe3
    
    U1 <- (A[1,1]*Ue1 + A[1,2]*Ue2 + A[1,3]*Ue3)  
    U2 <- (A[2,1]*Ue1 + A[2,2]*Ue2 + A[2,3]*Ue3)  
    U3 <- (A[3,1]*Ue1 + A[3,2]*Ue2 + A[3,3]*Ue3)  

    U1[U1<0] <- 0
    U2[U2<0] <- 0
    U3[U3<0] <- 0

    U1 <- sqrt(U1)
    U2 <- sqrt(U2)
    U3 <- sqrt(U3)

    return(cbind(-a1[1]*U1 - a2[1]*U2 - a3[1]*U3,
                 -a1[2]*U1 - a2[2]*U2 - a3[2]*U3,
                 -a1[3]*U1 - a2[3]*U2 - a3[3]*U3))
  }
}

#' Calibração de sensor 2D
#' @export
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


 

