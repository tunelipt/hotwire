\name{hotwire}
\alias{hotwire}
\title{Cálculo de velocidade e calibração de sondas de fio quente}
\description{Descrição geral das funções da biblioteca \code{hotwire}
  usadas no cálculo de velocidade e calibração de sondas de fio quente.}


\details{

  A biblioteca \code{hotwire} possui funções para cálculo de velocidade
  de sondas de fio quente. Aqui não existe nenhuma referência a como
  foram obtidos os dados do fio quente (tensões por exemplo).

  Por enquanto esta biblioteca limita-se a sondas unidimensionais
  (sondas simples). Para  sondas duplas ou triplas, os cálculos são mais
  complexos pois envolvem ângulos e outras contribuições. Este
  tratamento será implementado no futuro.

  A única função que não tem relação direta com anemometria de fio
  quente é a função \code{\link{v2temp}}. Esta função é a curva de
  calibração do termistor acoplado ao sistem 54N81 da DANTEC. Como este
  sensor é muito utilizando junto com o fio quente, resolveu-se colocar
  esta função nesta biblioteca.

  \strong{Procedimento de calibração}

  Agora é descrito como esta biblioteca ajuda na calibração de sondas de
  fio quente.

  Inicialmente é importante destacar que não se calibra uma sonda de fio
  quente e sim o conjunto sonda/eletrônica. A eletrônica pode ser mudada
  dependendo do caso (ver documentação do sistema de fio quente) e cada
  mudança requer uma nova calibração do conjunto.

  A calibração pode ser efetuada (no CMF/IPT) de três maneiras básicas:

  \enumerate{
    \item Sistema automático da DANTEC
    \item Sistema antigo da DANTEC
    \item Túnel de vento
  }

  O uso do sistema automático da DANTEC permite que velocidades sejam
  programadas e o valor pode ser lido diretamente. Não está claro como
  fazer a interface com o software da DANTEC mas uma leitura manual é
  possível.

  O sistema antigo da DANTEC permite fazer a calibração se estiver
  disponível um sistema de medição de pressões (scanivalve por
  exemplo). Com este sistema, pode-se estimar a velocidade no bocal e
  anotando estes valores esta biblioteca pode ser utilizada para fazer a
  calibração.

  Outra possibilidade é utilizar o túnel de vento pequeno para fazer a
  calibração. Novamente os dados devem ser anotados e depois as funções
  desta biblioteca R podem ser utilizadas para obter uma calibração.

  \strong{Uso de \code{hotwire} para obtenção de curvas de calibração}

  A função básica de calibração e a função \code{\link{hwCalibr}}. Dados
  os pontos de calibração e a configuração do sistema, esta função
  returna uma outra função que abstrai todos os aspectos da
  calibração. Esta outra função permite calcular a velocidade a partir
  da tensão medida e a temperatura ambiente na medição.

  Observe que a temperatura de calibração pode ser fornecida para cada
  ponto de calibração (correspondente a cada Ti e Vi) ou pode ser
  fornecido um valor único, válido para todos os pontos de calibração,
  se a temperatura não mudou muito durante a calibração.

  A temperatura de referência \code{T0} deve ser próxima da temperatura de uso mas
  pode ter qualquer valor na realidade.

  A temperatura do fio \code{Tw} é obtida a partir das características
  da sonda e a configuração da parte eletrônica. No caso do sistema
  54N81 da DANTEC, isto pode ser mudado utilizando os dipswitchs (DS) dentro
  do aparelho. A configuração dos DS vai depender da temperatura de
  operação desejada \code{Tw} e das resitências do fio, do suporte e
  cabos utilizados. Para configurar isso, favor procurar a documentação
  da DANTEC ou do sistema que está sendo utlizado. Para o sistem 54N81,
  a dantec fornece uma planilha excel que pode ser utilizada para
  configurar os DS. Ela pode ser encontrada em
  \url{http://www.dantecdynamics.com/Admin/Public/DWSDownload.aspx?File=%2fFiles%2ffiler%2fProducts+and+Services%2fGeneral+literature%2fMini-CTA_Excel_spreadsheet_v3A.XLS}.

    Uma cópia desta planilha pode ser encontrada no diretório retornado por:

    \code{.path.package('hotwire')}

    

  \strong{Cálculo de velocidade}

  Com o valor (função) retornado por \code{\link{hwCalibr}}, pode-se
  calcular diretamente a velocidade. Exemplo:

  \code{calibr <- hwCalibr(...)}
  \code{vel <- calibr(tensao, temperatura)}

  É lógico que esta biblioteca geralmente será utilizada em conjunto com
  algum sistema de aquisição de dados como por exemplo \pkg{Rnidaq}.

  \strong{O que falta ainda?}

  \enumerate{
    \item Sondas duplas e triplas (difícil)
    \item Banco de dados para armazenar dados de calibração (fácil)
  }
  
}

\seealso{
  \code{\link{v2temp}},
  \code{\link{fitKing}},
  \code{\link{kingVel}},
  \code{\link{tempCorr}},
  \code{\link{hwCalibr}},
  \code{\link{distrLog}}
}

\value{Função que calcula a a velocidade a partir da tensão}
