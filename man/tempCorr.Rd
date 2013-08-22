\name{tempCorr}
\alias{tempCorr}
\title{Correção de temperatura de sensor de fio quente.}
\usage{
  tempCorr(Ta, T0 = 25, Tw = 245)
}
\arguments{
  \item{Ta}{Temperatura ambiente.}

  \item{T0}{Temperatura de referência.}

  \item{Tw}{Temperatura de operação do fio quente.}
}
\value{
  Fator de correção que deve ser aplicada à tensão.
}
\description{
  Corrige a tensão lida em um sistema de fio quente como
  resultado da variação da temperatura ambiente em relação
  à temperatura ambiente de calibração.
}
\details{
  Caso a calibração seja feita a uma temperatura diferente
  da temperatura de operação, a tensão deve ser corrigida.
  Esta relação é baseada num artigo do Bearman na DISA
  Information no começo dos anos 70.
}

