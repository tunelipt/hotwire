\name{hotwire}
\alias{hotwire}
\title{C�lculo de velocidade e calibra��o de sondas de fio quente}
\description{Descri��o geral das fun��es da biblioteca hotwire
  usadas no c�lculo de velocidade e calibra��o de sondas de fio quente.}


\details{

  A biblioteca \pkg{hotwire} possui fun��es para c�lculo de velocidade
  de sondas de fio quente. Aqui n�o existe nenhuma refer�ncia a como
  foram obtidos os dados do fio quente (tens�es por exemplo).

  Por enquanto esta biblioteca limita-se a sondas unidimensionais
  (sondas simples). Para  sondas duplas ou triplas, os c�lculos s�o mais
  complexos pois envolvem �ngulos e outras contribui��es. Este
  tratamento ser� implementado no futuro.

  A �nica fun��o que n�o tem rela��o direta com anemometria de fio
  quente � a fun��o \code{\link{v2temp}}. Esta fun��o � a curva de
  calibra��o do termistor acoplado ao sistem 54N81 da DANTEC. Como este
  sensor � muito utilizando junto com o fio quente, resolveu-se colocar
  esta fun��o nesta biblioteca.

  \strong{Procedimento de calibra��o}

  Agora � descrito como esta biblioteca ajuda na calibra��o de sondas de
  fio quente.

  Inicialmente � importante destacar que n�o se calibra uma sonda de fio
  quente e sim o conjunto sonda/eletr�nica. A eletr�nica pode ser mudada
  dependendo do caso (ver documenta��o do sistema de fio quente) e cada
  mudan�a requer uma nova calibra��o do conjunto.

  A calibra��o pode ser efetuada (no CMF/IPT) de tr�s maneiras b�sicas:

  \enumerate{
    \item Sistema autom�tico da DANTEC
    \item Sistema antigo da DANTEC
    \item T�nel de vento
  }

  O uso do sistema autom�tico da DANTEC permite que velocidades sejam
  programadas e o valor pode ser lido diretamente. N�o est� claro como
  fazer a interface com o software da DANTEC mas uma leitura manual �
  poss�vel.

  O sistema antigo da DANTEC permite fazer a calibra��o se estiver
  dispon�vel um sistema de medi��o de press�es (scanivalve por
  exemplo). Com este sistema, pode-se estimar a velocidade no bocal e
  anotando estes valores esta biblioteca pode ser utilizada para fazer a
  calibra��o.

  Outra possibilidade � utilizar o t�nel de vento pequeno para fazer a
  calibra��o. Novamente os dados devem ser anotados e depois as fun��es
  desta biblioteca R podem ser utilizadas para obter uma calibra��o.

  \strong{Uso de \code{hotwire} para obten��o de curvas de calibra��o}

  A fun��o b�sica de calibra��o e a fun��o \code{\link{hwCalibr}}. Dados
  os pontos de calibra��o e a configura��o do sistema, esta fun��o
  returna uma outra fun��o que abstrai todos os aspectos da
  calibra��o. Esta outra fun��o permite calcular a velocidade a partir
  da tens�o medida e a temperatura ambiente na medi��o.

  Observe que a temperatura de calibra��o pode ser fornecida para cada
  ponto de calibra��o (correspondente a cada Ti e Vi) ou pode ser
  fornecido um valor �nico, v�lido para todos os pontos de calibra��o,
  se a temperatura n�o mudou muito durante a calibra��o.

  A temperatura de refer�ncia \code{T0} deve ser pr�xima da temperatura de uso mas
  pode ter qualquer valor na realidade.

  A temperatura do fio \code{Tw} � obtida a partir das caracter�sticas
  da sonda e a configura��o da parte eletr�nica. No caso do sistema
  54N81 da DANTEC, isto pode ser mudado utilizando os dipswitchs (DS) dentro
  do aparelho. A configura��o dos DS vai depender da temperatura de
  opera��o desejada \code{Tw} e das resit�ncias do fio, do suporte e
  cabos utilizados. Para configurar isso, favor procurar a documenta��o
  da DANTEC ou do sistema que est� sendo utlizado. Para o sistem 54N81,
  a dantec fornece uma planilha excel que pode ser utilizada para
  configurar os DS. Ela pode ser encontrada em no diret�rio retornado por:

  \code{.path.package('hotwire')}

    

  \strong{C�lculo de velocidade}

  Com o valor (fun��o) retornado por \code{\link{hwCalibr}}, pode-se
  calcular diretamente a velocidade. Exemplo:

  \code{calibr <- hwCalibr(...)}
  \code{vel <- calibr(tensao, temperatura)}

  � l�gico que esta biblioteca geralmente ser� utilizada em conjunto com
  algum sistema de aquisi��o de dados como por exemplo \pkg{Rnidaq}.

  \strong{O que falta ainda?}

  \enumerate{
    \item Sondas duplas e triplas (dif�cil)
    \item Banco de dados para armazenar dados de calibra��o (f�cil)
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

\value{Fun��o que calcula a a velocidade a partir da tens�o}
