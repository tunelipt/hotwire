# Implementa algumas opera�oes~comuns utilizadas
# em fio quente.


# O sistema da Dantec Multichannel CTA 54N81 possui um saida�para termistor.
# A curva de calibra�ao~obtida em 2005 e�a seguinte:
# Temp = -25.563*ln(tensao) + 4.984
# A cada ano ou dois anos seria interessante recalibrar o termistor e o
# circuito da DANTEC
v2temp <- function(V) -25.563 * log(V) + 4.984



# Caso a calibra��o seja feita a uma temperatura diferente da temperatura de
# opera��o, a tens�o deve ser corrigida. Esta rela��o � baseada num artigo
# do Bearman na DISA Information no come�o dos anos 70.
tempCorr <- function(Ta,T0=25, Tw=245){
	
  f <- sqrt( (Tw - T0) / (Tw - Ta) )
  return (f)
}


# Acha os coeficientes da lei de King generalizada para diferentes tens�es
# e velocidades. T^2 = A + B*V^n  (calcula os coeficientes A, B e n.
fitKing <- function(Ten, V){
  T2 <- Ten^2
  linear <- lm(T2 ~ I(V^0.5) )
  A <- linear[[1]][[1]]
  B <- linear[[1]][[2]]
  
  # Fazer o ajuste n�o linear
  erro <- function(p){
    sum( ( p[1] + p[2]*V^p[3] - T2)^2)
  }
  p <- c(A, B, 0.5)
  
  nonlin <- nlm(erro, p=p)
  
  x <- nonlin$estimate
  return(list(A=x[1], B=x[2], n=x[3]))
}

# Dados os pontos de calibra��o, cria uma fun��o que calcula a velocidade
# para diferentes tens�es.
kingVel <- function(Ten, V){
  fit <- fitKing(Ten, V)
  function(Ten) ((Ten^2-fit$A)/fit$B)^(1/fit$n)
}


  

# Cria uma sequ�ncia de n�meros com distribui��o log
distrLog <- function(xmin, xmax, N){

  lmin = log(xmin)
  lmax = log(xmax)
  
  lpos = seq(lmin, lmax, length=N)
  
  return(exp(lpos))
}



# Fun��o que realiza levanta a curva de calibra��o
hwCalibr <- function(tens, vel, temp, T0=25, Tw=245, curveFit=kingVel){
  
  f <- tempCorr(temp, T0, Tw)

  tens <- f * tens

  fitfun <- curveFit(tens, vel)

  # Criar a fun��o que aplica o ajuste de curva e calcula a velocide:
  function(tens, Ta=T0)
    fitfun(tens * tempCorr(Ta, T0, Tw))
}  
  