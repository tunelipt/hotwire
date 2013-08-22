\name{hwCalibr}
\alias{hwCalibr}
\title{Curva de calibração de sondas de fio quente}}
\usage{
  hwCalibr(tens, vel, temp, T0 = 25, Tw = 245,
    curveFit = kingVel)
}
\arguments{
  \item{tens}{Tensão medida na calibração}

  \item{Vel}{Velocidade de calibração}

  \item{temp}{Temperatura de calibração}

  \item{T0}{Temperatura de referência}

  \item{Tw}{Temperatura de operação do fio quente}

  \item{curveFit}{Função que faz o ajuste de curva}
}
\value{
  Função que calcula a a velocidade a partir da tensão
}
\description{
  Retorna função que calcula a velocidade medida no fio
  quente a partir da tensão e temperatura de medida.
}
\details{
  Esta função gerencia a calibração de uma sonda de fio
  quente. Na verdade o produto final desta função é uma
  função que calcula a velocidade medida pela sonda de fio
  quente / circuito eletrônico para diferentes tensões e
  temperaturas de operação.

  Ous seja, esta função retorna, a partir de parâmetros e
  medições de calibração uma função do tipo:

  \code{function(tens, Ta)}...

  Para fazer isso, a função primeiro corrige os valores de
  tensão para a condição de referência (para a mesma
  velocidade qual seria a tensão se a temperatura durante a
  medição fosse \code{T0}. Isto é feito chamando a função
  \code{\link{tempCorr}}.

  Com esta tensão nas condições de referência, chama-se a
  função \code{curveFit} para fazer um ajuste de curva.
  Este ajuste de curva é dado na forma de uma função que
  recebe a tensão e retorna a velocidade. O valor default
  deste parâmetro é \code{\link{kingVel}} que faz um ajuste
  utilizando a lei de king. Uma outra possibilidade é a
  função \code{\link{splinefun}} que interpola os dados
  utilizando uma spline cúbica. Qualquer função semelhante
  a \code{\link{splinefun}} pode ser utilizada.

  De posse desta função de interpolação/aproximação, é
  criada uma função que recebe dois parâmetros: a tensão e
  a temperatura de medição. Esta função primeiro corrige a
  tensão para as condições de referência (utilizando a
  função \code{\link{tempCorr}}) e calcula a velocidade
  utilizando o ajuste descrito na seção anterior.
}
\examples{
# Exemplo de velocidade e tensão
vel <- 1:15
tens <- sqrt(6 + 2*vel^0.43)

# Adicionar um pouco de ruído para simular dados experimentais:
velb <- vel + rnorm(length(vel), sd=0.3)

# Fazer o ajuste:
hw <- hwCalibr(tens, velb, temp=30, T0=25, Tw=245, curveFit=kingVel)

#Plotar os dados:
plot(tens, velb, xlab='Tensão (V)', ylab='Velocidade (m/s)',
main='Curva de Calibração de fio quente')

# Plotar o ajuste:

lines(tens, hw(tens, 30))
}
\seealso{
  \code{\link{fitKing}}, \code{\link{kingVel}},
  \code{\link{tempCorr}}, \code{\link{splinefun}}
}

