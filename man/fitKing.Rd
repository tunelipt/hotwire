\name{fitKing}
\alias{fitKing}
\title{Ajuste de curva segundo a lei de King generalizada}
\description{Acha os coeficientes da lei de King generalizada para uma
  sequência de medidas (tensão,velocidade)}

\usage{fit <- fitKing(Ten, Vel)}

\arguments{
  \item{Ten}{Tensão medida (pode ser corrigida)}
  \item{Vel}{Velocidade}
}

\details{
  Um dos modelos mais simples para o fio quente é a lei de king
  generalizada:
  \deqn{T^2 = A + B\cdot U^n}{T^2 = A+B.U^n}
  
  Esta é uma curva muitas vezes adequada para a representação de curvas
  de calibração. Observe que é importante que cada ponto obtido
  corresponda a uma mesma temperatura de referência. Se necessário é
  interessante utilizar a função \code{\link{tempCorr}} para se chegar a
  uma mesma temperatura de referência para todos os pontos.
  
}

\value{Lista com os coeficientes de ajuste A, B e n}

\examples{
# Exemplo de velocidade e tensão
vel <- 1:15
tens <- sqrt(6 + 2*vel^0.43)

# Adicionar um pouco de ruído para simular dados experimentais:
velb <- vel + rnorm(length(vel), sd=0.3)

# Fazer o ajuste:
fit <- fitKing(tens, velb)
cat("A = ", fit$A, '\nB = ', fit$B, '\nn = ', fit$n, '\n', sep='')

#Plotar os dados:
plot(tens, velb, xlab='Tensão (V)', ylab='Velocidade (m/s)',
main='Curva de Calibração de fio quente')

# Plotar o ajuste:

lines(tens, ((tens^2-fit$A)/fit$B)^(1/fit$n))
}


