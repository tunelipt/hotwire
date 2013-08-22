\name{distrLog}
\alias{distrLog}
\title{Sequencia de números com distribuição logarítimica}
\usage{
  distrLog(xmin, xmax, N)
}
\arguments{
  \item{xmin}{Ponto mínimo}

  \item{xmax}{Ponto máximo}

  \item{N}{Número total de números}
}
\description{
  Cria uma sequência de números que tem uma distribuição
  log. Pode ser utilizado para determinar as velocidades de
  calibração.
}
\examples{
A <- 6
B <- 2
n <- 0.45

u <- distrLog(1, 15, 10)

tens <- sqrt(B + A*u^n)

plot(u, tens, ty='b')
}

