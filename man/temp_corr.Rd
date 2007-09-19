\name{tempCorr}
\alias{tempCorr}
\title{Correção de temperatura de leitura de fio quente}
\description{Corrige a tensão lida num anemômetro de fio quente para a
  condição de calibração}

\usage{Vcorr <- tempCorr(Ta, T0, Tw}}

\arguments{
  \item{Ta}{Temperatura de medição (ambiente) em oC}
  \item{T0}{Temperatura de referência (temperatura de curva de
    calibração) em oC}
  \item{Tw}{Temperatura de operação do fio em oC}
}

\details{
  A medição de velocidade utilizando fio quente depende diretamente da
  transferência de calor entre o fio aquecido e o ambiente. Quando é
  feita uma calibração, geralmente a temperatura é aproximadamente
  constante, que chamaremos de T0.

  Quando o sensor é utilizado, a temperatura pode mudar
  consideravelmente, principalmente em túneis de vento onde a
  temperatura aumenta ao longo da operação. Quando a temperatura
  aumenta, para uma mesma velocidade de vento, a transferência de calor
  cai e portanto sem uma devida correção, a tensão lida é menor.

  Esta função realiza esta correção. Esta correção é baseada em um
  artigo de Bearman na DISA Information de 1971(?). A correção pode ser
  obtida a partir de correlações empíricas para o Nu e modelos
  qualitativos de funcionamento do fio quente. A correção obtida é:

  \deqn{V_{T_0} = V_T \sqrt{\frac{T_w - T_0}{T_w-T_a}}}{V_0 =
    V_T*((Tw-T0)/(Tw-Ta))^(1/2)}

  Esta função retorna a raíz quadrada da expressão anterior.
}

\value{Fator de correção da tensão}

\examples{
# Tensão lida no anemômetro
tens <- 3.56

f <- tempCorr(29.4, 25, 245)

tens_calibr <- f*tens
}


