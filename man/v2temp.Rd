\name{v2temp}
\alias{v2temp}
\title{Curva de calibra��o do termistor da dantec}
\description{Converte a tens�o lida no canal de temperatura do DANTEC
  54N81 na temperatura em oC}

\usage{temp <- v2temp(V) }

\arguments{
  \item{V}{Tens�o medida}
}

\details{
  Foi feita uma calibra��o do termistor e do circuito em conjunto em
  2005 utilizando os sistemas de calibra��o m�veis do laborat�rio de
  �leo e G�s do CMF/IPT.

  � recomendado que esta calibra��o seja feita a cada ano ou 2 anos.

  A curva obtida foi:

  \deqn{ T = -25.563 \ln V + 4.984}{T = -25.563*ln(V) + 4.984}
}

\value{Temperatura em oC}

\examples{
tens <- c(0.4, 0.43, 0.47, 0.5, 0.54, 0.6)
temperatura <- v2temp(tens)
plot(tens, temperatura, ty='b', xlab='Tens�o (V)',
ylab='Temperatura (oC)')
}






  
  

  