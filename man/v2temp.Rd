\name{v2temp}
\alias{v2temp}
\title{Cálculo de temperatura no Dantec Multichannel CTA 54N81}
\usage{
  v2temp(V)
}
\arguments{
  \item{V}{Tensão em Volts.}
}
\value{
  Temperatura em oC.
}
\description{
  O sistema da Dantec Multichannel CTA 54N81 possui um
  saida´para termistor. A curva de calibraçao~obtida em
  2005 e´a seguinte: Temp = -25.563*ln(tensao) + 4.984 A
  cada ano ou dois anos seria interessante recalibrar o
  termistor e o circuito da DANTEC
}

