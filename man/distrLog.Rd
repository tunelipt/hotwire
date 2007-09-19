\name{distrLog}
\alias{distrLog}
\title{Sequencia de n�meros com distribui��o logar�timica}
\description{Cria uma sequ�ncia de n�meros que tem uma distribui��o
  log. Pode ser utilizado para determinar as velocidades de calibra��o.}

\usage{x <- distrLog(xmin, xmax, N)}

\arguments{
  \item{xmin}{Ponto m�nimo}
  \item{xmax}{Ponto m�ximo}
  \item{N}{N�mero total de n�meros}
}


\value{Vetor de n�meros}

\examples{
A <- 6
B <- 2
n <- 0.45

u <- distrLog(1, 15, 10)

tens <- sqrt(B + A*u^n)

plot(u, tens, ty='b')
}

