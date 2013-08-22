\name{kingVel}
\alias{kingVel}
\title{Ajuste de curva segundo a lei de King generalizada}
\usage{
  kingVel(Ten, V)
}
\arguments{
  \item{Ten}{Tensão de calibração.}

  \item{V}{Velocidade de calibração.}
}
\value{
  Função que estima a velocidade a partir da tensão.
}
\description{
  Retorna função que calcula a velocidade a partir da
  tensão
}
\details{
  Um dos modelos mais simples para o fio quente é a lei de
  king generalizada: \deqn{T^2 = A + B\cdot U^n}{T^2 =
  A+B.U^n}

  Esta é uma curva muitas vezes adequada para a
  representação de curvas de calibração. Observe que é
  importante que cada ponto obtido corresponda a uma mesma
  temperatura de referência. Se necessário é interessante
  utilizar a função \code{\link{tempCorr}} para se chegar a
  uma mesma temperatura de referência para todos os pontos.

  Esta função chama a função \code{\link{fitKing}} para
  calcular os coeficientes de ajuste. Com estes
  coeficientes, esta função cria uma função que calcula a
  velocidade a partir da tensão medida e corrigida para a
  referência.
}
\examples{
# Exemplo de velocidade e tensão
vel <- 1:15
tens <- sqrt(6 + 2*vel^0.43)

# Adicionar um pouco de ruído para simular dados experimentais:
velb <- vel + rnorm(length(vel), sd=0.3)

# Fazer o ajuste:
fitfun <- kingVel(tens, velb)

#Plotar os dados:
plot(tens, velb, xlab='Tensão (V)', ylab='Velocidade (m/s)',
main='Curva de Calibração de fio quente')

# Plotar o ajuste:

lines(tens, fitfun(tens))
}
\seealso{
  \code{\link{kingVel}}
}

