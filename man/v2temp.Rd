\name{v2temp}
\alias{v2temp}
\title{Curva de calibração do termistor da dantec}
\description{Converte a tensão lida no canal de temperatura do DANTEC
  54N81 na temperatura em oC}

\usage{temp <- v2temp(V) }

\arguments{
  \item{V}{Tensão medida}
}

\details{
  Foi feita uma calibração do termistor e do circuito em conjunto em
  2005 utilizando os sistemas de calibração móveis do laboratório de
  Óleo e Gás do CMF/IPT.

  É recomendado que esta calibração seja feita a cada ano ou 2 anos.

  A curva obtida foi:

  \deqn{ T = -25.563 \ln V + 4.984}{T = -25.563*ln(V) + 4.984}
}

\value{Temperatura em oC}

\examples{
tens <- c(0.4, 0.43, 0.47, 0.5, 0.54, 0.6)
temperatura <- v2temp(tens)
plot(tens, temperatura, ty='b', xlab='Tensão (V)',
ylab='Temperatura (oC)')
}






  
  

  