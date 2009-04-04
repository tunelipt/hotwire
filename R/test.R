test2D <- function(){

  support <- hwSupport(contacts=2)
  cable <- hwCable()
  signal <- hwSignalConditioner()
  bridge <- hwBridgeConfig()+
  w1 <- hwWire(3.3, X=pi/4, Y=-pi/4)
  w2 <- hwWire(3.4, X=3*pi/4, Y=pi/4)

  probe <- hwProbe(support, cable, bridge, signal, w1, w2)
  U <- 1:20

  E1 <- sqrt(2 + 1 * U^0.5)
  E2 <- sqrt(2.2 + 1.1*U^0.5)

  calibr <- hwCalibrProbe2D(probe, E1, U, E2, U)
  return(calibr)
}



test3D <- function(){
  d2r <- function(d) d*pi/180

  support <- hwSupport(contacts=2)
  cable <- hwCable()
  signal <- hwSignalConditioner()
  bridge <- hwBridgeConfig()
  
  w1 <- hwWire(R=6,   k=0.04, h=1.2, X=d2r(125.264), Y=d2r(45),  Z=d2r(114.094))
  w2 <- hwWire(R=6.1, k=0.04, h=1.2, X=d2r(125.264), Y=d2r(135), Z=d2r(114.094))
  w3 <- hwWire(R=6.2, k=0.04, h=1.2, X=d2r(125.264), Y=d2r(90),  Z=d2r(35.264))
  probe <- hwProbe(support, cable, bridge, signal, w1, w2, w3)
  U <- 1:20

  E1 <- sqrt(2 + 1 * U^0.5)
  E2 <- sqrt(2.2 + 1.1*U^0.5)
  E3 <- sqrt(2.1 + 1.2*U^0.5)

  calibr <- hwCalibrProbe3D(probe, E1, U, E2, U, E3, U)
  return(calibr)
}
