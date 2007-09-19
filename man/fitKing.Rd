\name{fitKing}
\alias{fitKing}
\title{Ajuste de curva segundo a lei de King generalizada}
\description{Acha os coeficientes da lei de King generalizada para uma
  sequ�ncia de medidas (tens�o,velocidade)}

\usage{fit <- fitKing(Ten, Vel)}

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
  
}

\value{Lista com os coeficientes de ajuste A, B e n}

\examples{
# Exemplo de velocidade e tens�o
vel <- 1:15
tens <- sqrt(6 + 2*vel^0.43)

# Adicionar um pouco de ru�do para simular dados experimentais:
velb <- vel + rnorm(length(vel), sd=0.3)

# Fazer o ajuste:
fit <- fitKing(tens, velb)
cat("A = ", fit$A, '\nB = ', fit$B, '\nn = ', fit$n, '\n', sep='')

#Plotar os dados:
plot(tens, velb, xlab='Tens�o (V)', ylab='Velocidade (m/s)',
main='Curva de Calibra��o de fio quente')

# Plotar o ajuste:

lines(tens, ((tens^2-fit$A)/fit$B)^(1/fit$n))
}

