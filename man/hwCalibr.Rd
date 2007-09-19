\name{hwCalibr}
\alias{hwCalibr}
\title{Curva de calibra��o de sondas de fio quente}
\description{Retorna fun��o que calcula a velocidade medida no fio
  quente a partir da tens�o e temperatura de medida}

\usage{calibr <- hwCalibr(tens, vel, temp, T0, Tw, curveFit)}

\arguments{
  \item{tens}{Tens�o medida na calibra��o}
  \item{Vel}{Velocidade de calibra��o}
  \item{temp}{Temperatura de calibra��o}
  \item{T0}{Temperatura de refer�ncia}
  \item{Tw}{Temperatura de opera��o do fio quente}
  \item{curveFit}{Fun��o que faz o ajuste de curva}
}

\details{

  Esta fun��o gerencia a calibra��o de uma sonda de fio quente. Na
  verdade o produto final desta fun��o � uma fun��o que calcula a
  velocidade medida pela sonda de fio quente / circuito eletr�nico para
  diferentes tens�es e temperaturas de opera��o.

  Ous seja, esta fun��o retorna, a partir de par�metros e medi��es de
  calibra��o uma fun��o do tipo:

  \code{function(tens, Ta)}...

  Para fazer isso, a fun��o primeiro corrige os valores de tens�o para a
  condi��o de refer�ncia (para a mesma velocidade qual seria a tens�o se
  a temperatura durante a medi��o fosse \code{T0}. Isto � feito chamando
  a fun��o \code{\link{tempCorr}}.

  Com esta tens�o nas condi��es de refer�ncia, chama-se a fun��o
  \code{curveFit} para fazer um ajuste de curva. Este ajuste de curva �
  dado na forma de uma fun��o que recebe a tens�o e retorna a
  velocidade. O valor default deste par�metro � \code{\link{kingVel}}
  que faz um ajuste utilizando a lei de king. Uma outra possibilidade �
  a fun��o \code{\link{splinefun}} que interpola os dados utilizando uma
  spline c�bica. Qualquer fun��o semelhante a \code{\link{splinefun}}
  pode ser utilizada.

  De posse desta fun��o de interpola��o/aproxima��o, � criada uma fun��o
  que recebe dois par�metros: a tens�o e a temperatura de medi��o. Esta
  fun��o primeiro corrige a tens�o para as condi��es de refer�ncia
  (utilizando a fun��o \code{\link{tempCorr}}) e calcula a velocidade
  utilizando o ajuste descrito na se��o anterior.
  

  
}

\seealso{
  \code{\link{fitKing}}, \code{\link{kingVel}}, \code{\link{tempCorr}}, \code{\link{splinefun}}
}

\value{Fun��o que calcula a a velocidade a partir da tens�o}

\examples{
# Exemplo de velocidade e tens�o
vel <- 1:15
tens <- sqrt(6 + 2*vel^0.43)

# Adicionar um pouco de ru�do para simular dados experimentais:
velb <- vel + rnorm(length(vel), sd=0.3)

# Fazer o ajuste:
hw <- hwCalibr(tens, velb, temp=30, T0=25, Tw=245, curveFit=kingVel)

#Plotar os dados:
plot(tens, velb, xlab='Tens�o (V)', ylab='Velocidade (m/s)',
main='Curva de Calibra��o de fio quente')

# Plotar o ajuste:

lines(tens, hw(tens, 30))
}

