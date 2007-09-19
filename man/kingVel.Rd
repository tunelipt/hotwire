\name{kingVel}
\alias{kingVel}
\title{Ajuste de curva segundo a lei de King generalizada}
\description{Retorna fun��o que calcula a velocidade a partir da tens�o}

\usage{fitfun <- kingVel(Ten, Vel)}

\arguments{
  \item{Ten}{Tens�o medida (pode ser corrigida)}
  \item{Vel}{Velocidade}
}

\details{
  Um dos modelos mais simples para o fio quente � a lei de king
  generalizada:
  \deqn{T^2 = A + B\cdot U^n}{T^2 = A+B.U^n}
  
  Esta � uma curva muitas vezes adequada para a representa��o de curvas
  de calibra��o. Observe que � importante que cada ponto obtido
  corresponda a uma mesma temperatura de refer�ncia. Se necess�rio �
  interessante utilizar a fun��o \code{\link{tempCorr}} para se chegar a
  uma mesma temperatura de refer�ncia para todos os pontos.

  Esta fun��o chama a fun��o \code{\link{fitKing}} para calcular os
  coeficientes de ajuste. Com estes coeficientes, esta fun��o cria uma
  fun��o que calcula a velocidade a partir da tens�o medida e corrigida
  para a refer�ncia.
  
}

\seealso{
  \code{\link{fitKing}}
}

\value{Fun��o que calcula a a velocidade a partir da tens�o}

\examples{
# Exemplo de velocidade e tens�o
vel <- 1:15
tens <- sqrt(6 + 2*vel^0.43)

# Adicionar um pouco de ru�do para simular dados experimentais:
velb <- vel + rnorm(length(vel), sd=0.3)

# Fazer o ajuste:
fitfun <- kingVel(tens, velb)

#Plotar os dados:
plot(tens, velb, xlab='Tens�o (V)', ylab='Velocidade (m/s)',
main='Curva de Calibra��o de fio quente')

# Plotar o ajuste:

lines(tens, fitfun(tens))
}

