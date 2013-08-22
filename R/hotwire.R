# Implementa algumas operaçoes~comuns utilizadas
# em fio quente.

#' Cálculo de temperatura no Dantec Multichannel CTA 54N81
#'
#' O sistema da Dantec Multichannel CTA 54N81 possui um saida´para termistor.
#' A curva de calibraçao~obtida em 2005 e´a seguinte:
#' Temp = -25.563*ln(tensao) + 4.984
#' A cada ano ou dois anos seria interessante recalibrar o termistor e o
#' circuito da DANTEC
#'
#' @param V Tensão em Volts.
#' @return Temperatura em oC.
#' @export
v2temp <- function(V) -25.563 * log(V) + 4.984


#' Correção de temperatura de sensor de fio quente.
#'
#' Corrige a tensão lida em um sistema de fio quente
#' como resultado da variação da temperatura ambiente
#' em relação à temperatura ambiente de calibração.
#'
#' Caso a calibração seja feita a uma temperatura diferente da temperatura de
#' operação, a tensão deve ser corrigida. Esta relação é baseada num artigo
#' do Bearman na DISA Information no começo dos anos 70.
#'
#' @param Ta Temperatura ambiente.
#' @param T0 Temperatura de referência.
#' @param Tw Temperatura de operação do fio quente.
#' @return Fator de correção que deve ser aplicada à tensão.
#' @export
tempCorr <- function(Ta,T0=25, Tw=245){
	
  f <- sqrt( (Tw - T0) / (Tw - Ta) )
  return (f)
}


#' Ajuste de curva segundo a lei de King generalizada
#'
#' Acha os coeficientes da lei de King generalizada para uma
#' sequência de medidas (tensão,velocidade)
#'
#'  Um dos modelos mais simples para o fio quente é a lei de king
#'  generalizada:
#'  \deqn{T^2 = A + B\cdot U^n}{T^2 = A+B.U^n}
#'  
#'  Esta é uma curva muitas vezes adequada para a representação de curvas
#'  de calibração. Observe que é importante que cada ponto obtido
#'  corresponda a uma mesma temperatura de referência. Se necessário é
#'  interessante utilizar a função \code{\link{tempCorr}} para se chegar a
#'  uma mesma temperatura de referência para todos os pontos.
#'
#' @param Ten Tensão de calibração.
#' @param V Velocidade de calibração.
#' @return Lista com os parâmetros da Lei de King generalizada A, B e n.
#' @seealso \code{\link{kingVel}}
#' @examples
#' # Exemplo de velocidade e tensão
#' vel <- 1:15
#' tens <- sqrt(6 + 2*vel^0.43)
#' 
#' # Adicionar um pouco de ruído para simular dados experimentais:
#' velb <- vel + rnorm(length(vel), sd=0.3)
#' 
#' # Fazer o ajuste:
#' fit <- fitKing(tens, velb)
#' cat("A = ", fit$A, '\nB = ', fit$B, '\nn = ', fit$n, '\n', sep='')
#' 
#' #Plotar os dados:
#' plot(tens, velb, xlab='Tensão (V)', ylab='Velocidade (m/s)',
#' main='Curva de Calibração de fio quente')
#' 
#' # Plotar o ajuste:
#' 
#' lines(tens, ((tens^2-fit$A)/fit$B)^(1/fit$n))
#' @export
fitKing <- function(Ten, V){
  T2 <- Ten^2
  linear <- lm(T2 ~ I(V^0.5) )
  A <- linear[[1]][[1]]
  B <- linear[[1]][[2]]
  
  # Fazer o ajuste não linear
  erro <- function(p){
    sum( ( p[1] + p[2]*V^p[3] - T2)^2)
  }
  p <- c(A, B, 0.5)
  
  nonlin <- nlm(erro, p=p)
  
  x <- nonlin$estimate
  return(list(A=x[1], B=x[2], n=x[3]))
}

#' Ajuste de curva segundo a lei de King generalizada
#'
#' Retorna função que calcula a velocidade a partir da tensão
#'
#'  Um dos modelos mais simples para o fio quente é a lei de king
#'  generalizada:
#'  \deqn{T^2 = A + B\cdot U^n}{T^2 = A+B.U^n}
#'  
#'  Esta é uma curva muitas vezes adequada para a representação de curvas
#'  de calibração. Observe que é importante que cada ponto obtido
#'  corresponda a uma mesma temperatura de referência. Se necessário é
#'  interessante utilizar a função \code{\link{tempCorr}} para se chegar a
#'  uma mesma temperatura de referência para todos os pontos.
#'
#'  Esta função chama a função \code{\link{fitKing}} para calcular os
#'  coeficientes de ajuste. Com estes coeficientes, esta função cria uma
#'  função que calcula a velocidade a partir da tensão medida e corrigida
#'  para a referência.
#' @param Ten Tensão de calibração.
#' @param V Velocidade de calibração.
#' @return Função que estima a velocidade a partir da tensão.
#' @seealso \code{\link{kingVel}}
#' @examples
#' # Exemplo de velocidade e tensão
#' vel <- 1:15
#' tens <- sqrt(6 + 2*vel^0.43)
#' 
#' # Adicionar um pouco de ruído para simular dados experimentais:
#' velb <- vel + rnorm(length(vel), sd=0.3)
#' 
#' # Fazer o ajuste:
#' fitfun <- kingVel(tens, velb)
#' 
#' #Plotar os dados:
#' plot(tens, velb, xlab='Tensão (V)', ylab='Velocidade (m/s)',
#' main='Curva de Calibração de fio quente')
#' 
#' # Plotar o ajuste:
#' 
#' lines(tens, fitfun(tens))
#' @export
kingVel <- function(Ten, V){
  fit <- fitKing(Ten, V)
  function(Ten) ((Ten^2-fit$A)/fit$B)^(1/fit$n)
}


#' Ajuste polinomial.
#'
#' Cria função  que usa o método de mínimos quadrados para fazer um ajuste polinomial
#'
#' @param degree Grau do polinômio de ajuste.
#' @return Função que retorna returna uma função de ajuste polinomial.
#' @export
makePolyFit <- function(degree){
  function(x,y){
    fit <- lm(y ~ poly(x, degree))
    function(x)
      as.vector(predict.lm(fit, newdata=data.frame(x)))
  }
}


#' Sequencia de números com distribuição logarítimica
#'
#' Cria uma sequência de números que tem uma distribuição
#' log. Pode ser utilizado para determinar as velocidades de calibração.
#' @param xmin Ponto mínimo
#' @param xmax Ponto máximo
#' @param N Número total de números
#' @examples
#' A <- 6
#' B <- 2
#' n <- 0.45
#' 
#' u <- distrLog(1, 15, 10)
#' 
#' tens <- sqrt(B + A*u^n)
#' 
#' plot(u, tens, ty='b')
#' @export
distrLog <- function(xmin, xmax, N){

  lmin = log(xmin)
  lmax = log(xmax)
  
  lpos = seq(lmin, lmax, length=N)
  
  return(exp(lpos))
}



#' Curva de calibração de sondas de fio quente
#'
#' Retorna função que calcula a velocidade medida no fio
#' quente a partir da tensão e temperatura de medida.
#'
#' 
#'   Esta função gerencia a calibração de uma sonda de fio quente. Na
#'   verdade o produto final desta função é uma função que calcula a
#'   velocidade medida pela sonda de fio quente / circuito eletrônico para
#'   diferentes tensões e temperaturas de operação.
#' 
#'   Ous seja, esta função retorna, a partir de parâmetros e medições de
#'   calibração uma função do tipo:
#' 
#'   \code{function(tens, Ta)} ...
#' 
#'   Para fazer isso, a função primeiro corrige os valores de tensão para a
#'   condição de referência (para a mesma velocidade qual seria a tensão se
#'   a temperatura durante a medição fosse \code{T0}. Isto é feito chamando
#'   a função \code{\link{tempCorr}}.
#' 
#'   Com esta tensão nas condições de referência, chama-se a função
#'   \code{curveFit} para fazer um ajuste de curva. Este ajuste de curva é
#'   dado na forma de uma função que recebe a tensão e retorna a
#'   velocidade. O valor default deste parâmetro é \code{\link{kingVel}}
#'   que faz um ajuste utilizando a lei de king. Uma outra possibilidade é
#'   a função \code{\link{splinefun}} que interpola os dados utilizando uma
#'   spline cúbica. Qualquer função semelhante a \code{\link{splinefun}}
#'   pode ser utilizada.
#' 
#'   De posse desta função de interpolação/aproximação, é criada uma função
#'   que recebe dois parâmetros: a tensão e a temperatura de medição. Esta
#'   função primeiro corrige a tensão para as condições de referência
#'   (utilizando a função \code{\link{tempCorr}}) e calcula a velocidade
#'   utilizando o ajuste descrito na seção anterior.
#'
#'  @param tens Tensão medida na calibração
#'  @param Vel Velocidade de calibração
#'  @param temp Temperatura de calibração
#'  @param T0 Temperatura de referência
#'  @param Tw Temperatura de operação do fio quente
#'  @param curveFit Função que faz o ajuste de curva
#'  @return Função que calcula a a velocidade a partir da tensão
#' @seealso \code{\link{fitKing}}, \code{\link{kingVel}}, \code{\link{tempCorr}}, \code{\link{splinefun}}
#' @examples
#' # Exemplo de velocidade e tensão
#' vel <- 1:15
#' tens <- sqrt(6 + 2*vel^0.43)
#' 
#' # Adicionar um pouco de ruído para simular dados experimentais:
#' velb <- vel + rnorm(length(vel), sd=0.3)
#' 
#' # Fazer o ajuste:
#' hw <- hwCalibr(tens, velb, temp=30, T0=25, Tw=245, curveFit=kingVel)
#' 
#' #Plotar os dados:
#' plot(tens, velb, xlab='Tensão (V)', ylab='Velocidade (m/s)',
#' main='Curva de Calibração de fio quente')
#' 
#' # Plotar o ajuste:
#' 
#' lines(tens, hw(tens, 30))
#' @export
hwCalibr <- function(tens, vel, temp, T0=25, Tw=245, curveFit=kingVel){
  
  f <- tempCorr(temp, T0, Tw)

  tens <- f * tens

  fitfun <- curveFit(tens, vel)

  # Criar a função que aplica o ajuste de curva e calcula a velocide:
  function(tens, Ta=T0)
    fitfun(tens * tempCorr(Ta, T0, Tw))
}  
  
