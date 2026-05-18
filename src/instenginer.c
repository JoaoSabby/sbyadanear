/*
 * Projeto: adanear
 * Autor: Joao Batista Goncalves de Brito
 * 
 * Justificativa Arquitetural para Implementacao em C:
 * Esta biblioteca nativa foi desenvolvida para transpor gargalos criticos de 
 * performance inerentes a execucao de loops massivos no ambiente R puro. 
 * Algoritmos de superamostragem como o ADASYN e rotinas de padronizacao 
 * exigem manipulacao intensiva de tensores multidimensionais. 
 * 
 * A integracao em baixo nivel fornece controle absoluto sobre a alocacao da heap, 
 * permitindo operacoes diretas em memoria via ponteiros escalares sem acionar 
 * copias intermediarias de objetos. Alem de garantir tempos de execucao drasticamente 
 * menores em pipelines de modelagem preditiva, a base em C blinda a estabilidade 
 * matematica ao utilizar variaveis de precisao estendida no tratamento de 
 * variancias proximas a zero, isolando as operacoes do coletor de lixo da maquina virtual R.
 *
 * Referencias Historicas Fundamentais:
 * 1. He, H., Bai, Y., Garcia, E. A., & Li, S. (2008). ADASYN: Adaptive synthetic 
 *    sampling approach for imbalanced learning. IEEE International Joint Conference 
 *    on Neural Networks. (Base fundacional do algoritmo de superamostragem implementado).
 *
 * 2. Knuth, D. E. (1997). The Art of Computer Programming, Volume 2: Seminumerical 
 *    Algorithms. Addison Wesley. (Base fundacional historica para a estabilidade 
 *    numerica e algoritmos seguros de calculo de variancia em ponto flutuante).
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <R_ext/Visibility.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <math.h>
#include <limits.h>
#include <float.h>

/* FCONE: macro de R-API para passar comprimentos de strings Fortran nas
 * chamadas a BLAS. Definida em Rconfig.h em R >= 3.6; fornece um fallback
 * vazio em versoes mais antigas para manter portabilidade. */
#ifndef FCONE
# define FCONE
#endif

/**
 * @brief Invoca verificacao de interrupcao do usuario no ambiente R.
 *
 * Previne travamentos durante execucoes prolongadas na camada C, repassando 
 * o sinal para a API do R avaliar quebras manuais via teclado.
 *
 * @return SEXP Retorna R_NilValue incondicionalmente.
 */
SEXP OU_CheckUserInterruptC(void){
  
  /* Repassa o sinal para a API do R avaliar quebras manuais via teclado */
  R_CheckUserInterrupt();
  return R_NilValue;
}

/**
 * @brief Valida se o objeto SEXP fornecido e uma matriz de precisao dupla (REALSXP).
 *
 * Interrompe a execucao com erro caso o tipo ou formato divirjam do tensor 2D de reais.
 *
 * @param x O objeto SEXP a ser inspecionado.
 * @param name O nome da variavel injetado na formatacao do log de erro.
 */
static void require_real_matrix(SEXP x, const char *name){
  
  /* Trava a execucao caso a entrada divirja do formato tensor 2D de reais */
  if (!isReal(x) || !isMatrix(x)) {
    error("O parametro '%s' deve ser uma matrix double", name);
  }
}

/**
 * @brief Calcula media e desvio padrao amostral por coluna para normalizacao.
 *
 * Implementa o algoritmo de duas passagens com tipos long double para
 * evitar cancelamento catastrofico no calculo da variancia.
 *
 * @param xMatrix SEXP (matriz numerica double) contendo a base de dados.
 * @return SEXP Lista R contendo vetores numericos "centers" e "scales".
 */
SEXP OU_ComputeZScoreParamsC(SEXP xMatrix){
  
  /* Valida coerencia da estrutura de dados de entrada */
  require_real_matrix(xMatrix, "xMatrix");
  
  /* Extrai metadados do objeto para capturar quantidades de linhas e colunas.
   * Mantemos a captura como int porque R_DimSymbol entrega INTSXP, mas
   * promovemos para R_xlen_t antes de qualquer multiplicacao de offset para
   * evitar overflow em matrizes maiores que INT_MAX em uma das dimensoes. */
  SEXP dims = getAttrib(xMatrix, R_DimSymbol);
  const R_xlen_t n = (R_xlen_t) INTEGER(dims)[0];
  const R_xlen_t p = (R_xlen_t) INTEGER(dims)[1];

  /* Impede divisao por zero no calculo da variancia amostral n menos 1 */
  if(n < 2){
    error("'xMatrix' deve conter ao menos duas linhas");
  }
  
  /* Reserva blocos de memoria para media e desvio e inibe acao do Garbage Collector do R */
  SEXP centers = PROTECT(allocVector(REALSXP, p));
  SEXP scales = PROTECT(allocVector(REALSXP, p));
  
  /* Instancia ponteiros C brutos apontando para os enderecos mapeados das estruturas R */
  const double *x = REAL(xMatrix);
  double *mu = REAL(centers);
  double *sd = REAL(scales);
  
  /* Itera pelas colunas mapeando indices contiguos de memoria column-major */
  for(R_xlen_t j = 0; j < p; ++j){
    /* Previne sobrecarga no processador verificando interrupcoes apenas a cada 16 ciclos */
    if((j & 15) == 0){
      R_CheckUserInterrupt();
    }

    /* Posiciona o ponteiro de leitura no byte zero da coluna atual */
    const double *col = x + j * n;
    long double sum = 0.0L;

    /* Executa primeira passagem agregando valores para derivar o momento inicial absoluto */
    for(R_xlen_t i = 0; i < n; ++i){
      sum += (long double) col[i];
    }

    /* Define o parametro de centralizacao e grava no vetor protegido */
    const long double mean = sum / (long double) n;
    mu[j] = (double) mean;

    long double ss = 0.0L;

    /* Executa segunda passagem acumulando quadrados residuais estritamente sobre residuos centrados */
    for (R_xlen_t i = 0; i < n; ++i) {
      const long double centered = (long double) col[i] - mean;
      ss += centered * centered;
    }

    /* Conclui calculo da raiz do espalhamento da amostra e grava no vetor de escalas */
    sd[j] = sqrt((double) (ss / (long double) (n - 1)));
  }
  
  /* Constroi dicionario de nomes chaves para empacotamento na interface R */
  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("centers"));
  SET_STRING_ELT(names, 1, mkChar("scales"));
  
  /* Define bloco raiz tipo lista agrupando os resultados dos dois parametros */
  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, centers);
  SET_VECTOR_ELT(out, 1, scales);
  
  /* Vincula os rotulos string aos indices da lista recem criada */
  setAttrib(out, R_NamesSymbol, names);
    
  /* Decrementa contador de protecao em pilha liberando os 4 blocos locais gerados */
  UNPROTECT(4);
  return out;
}

/**
 * @brief Executa a transformacao Z-Score ou reversao dimensional em matriz double.
 *
 * Utiliza logica de ramificacao externa aos lacos internos para ganho de performance,
 * reposicionando os dados ao redor da media especificada e ajustando pelo desvio.
 *
 * @param xMatrix SEXP (matriz numerica double) com os dados alvo.
 * @param centers SEXP (vetor numerico double) contendo o momento inicial absoluto (media).
 * @param scales SEXP (vetor numerico double) contendo a raiz do espalhamento (desvio).
 * @param reverse SEXP (booleano logico). Valor TRUE restaura a escala nativa.
 * @return SEXP Matriz numerica R contendo a base transformada.
 */
SEXP OU_ApplyZScoreC(SEXP xMatrix, SEXP centers, SEXP scales, SEXP reverse){
  
  /* Valida existencia de matriz alvo e arrays de normalizacao pre ajustados */
  require_real_matrix(xMatrix, "xMatrix");
  
  if(!isReal(centers) || !isReal(scales)){
    error("'centers' e 'scales' devem ser vetores double");
  }
  
  /* Carrega metrica dimensional da base para validacao cruzada */
  SEXP dims = getAttrib(xMatrix, R_DimSymbol);
  const R_xlen_t n = (R_xlen_t) INTEGER(dims)[0];
  const R_xlen_t p = (R_xlen_t) INTEGER(dims)[1];

  /* Assegura conformidade de espaco vetorial de colunas x vetores de distribuicao */
  if(XLENGTH(centers) != p || XLENGTH(scales) != p){
    error("'centers' e 'scales' devem ter comprimento igual ao numero de colunas");
  }
  
  /* Cria espelho da matriz alvo alocando novo bloco de memoria homogeneo. */
  if(n > INT_MAX || p > INT_MAX){
    error("'xMatrix' excede o limite de dimensoes suportado por allocMatrix");
  }
  SEXP out = PROTECT(allocMatrix(REALSXP, (int) n, (int) p));
  
  /* Estabelece vinculacoes de enderecos para varredura C padrao */
  const double *x = REAL(xMatrix);
  const double *mu = REAL(centers);
  const double *sd = REAL(scales);
  double *y = REAL(out);
  
  /* Traduz sinalizador booleano logico do R em inteiro avaliavel em C */
  const int do_reverse = asLogical(reverse);
  
  for(R_xlen_t j = 0; j < p; ++j){
    /* Restringe consultas de hook de sistema operacional */
    if((j & 15) == 0){
      R_CheckUserInterrupt();
    }

    /* Pre calcula salto contiguo de bytes para a base do eixo y */
    const R_xlen_t offset = j * n;
    const double center = mu[j];
    const double scale = sd[j];

    /* Defesa em profundidade: sby_validate_scaling_info ja garante scale > 0,
     * mas duplicamos a checagem aqui porque OU_ApplyZScoreC pode ser invocado
     * com scaling_info fornecido externamente (sby_input_already_scaled). */
    if(!do_reverse && (!R_FINITE(scale) || scale <= 0.0)){
      error("'scales[%lld]' deve ser positivo e finito", (long long) (j + 1));
    }
    if(!R_FINITE(center)){
      error("'centers[%lld]' deve ser finito", (long long) (j + 1));
    }

    if(do_reverse){
      /* Restaura posicao real baseando se em distribuicao escalada pre formatada */
      for(R_xlen_t i = 0; i < n; ++i){
        y[offset + i] = x[offset + i] * scale + center;
      }
    } else {
      /* Normaliza vetor realocando ao redor de zero em desvios unificados unitarios */
      for(R_xlen_t i = 0; i < n; ++i){
        y[offset + i] = (x[offset + i] - center) / scale;
      }
    }
  }
  
  /* Copia rotulos e nomes contidos nas dimensoes de linha e coluna alvo */
  setAttrib(out, R_DimNamesSymbol, getAttrib(xMatrix, R_DimNamesSymbol));
  
  /* Retira ultimo encapsulamento pendente antes de despachar o resultado */
  UNPROTECT(1);
  return out;
}

/**
 * @brief Gera dados sinteticos aplicando interpolacao topologica (ADASYN).
 *
 * Expande a representatividade da classe alvo minoritaria atraves da criacao 
 * estocastica de novas ocorrencias posicionadas no plano cartesiano entre 
 * a observacao base e seus k-vizinhos diretos.
 *
 * @param minorityMatrix SEXP (matriz double) Base de features da classe limitante.
 * @param minorityNeighborIndex SEXP (matriz integer) Indices referenciais da malha vizinha.
 * @param syntheticPerRow SEXP (vetor integer) Carga de reproducao individualizada por linha.
 * @return SEXP Matriz double compreendendo unicamente a nova massa de dados gerada.
 */
SEXP OU_GenerateSyntheticAdasynC(SEXP minorityMatrix, SEXP minorityNeighborIndex, SEXP syntheticPerRow){
  /* Valida restricoes primarias de construcao da rede relacional vizinha */
  require_real_matrix(minorityMatrix, "minorityMatrix");
  
  if(!isInteger(minorityNeighborIndex) || !isMatrix(minorityNeighborIndex)){
    error("'minorityNeighborIndex' deve ser uma matrix integer");
  }
  
  if(!isInteger(syntheticPerRow)){
    error("'syntheticPerRow' deve ser um vetor integer");
  }
  
  /* Captura e indexa dimensoes de todos os tensores envolvidos na multiplicacao */
  SEXP minorityDims = getAttrib(minorityMatrix, R_DimSymbol);
  SEXP neighborDims = getAttrib(minorityNeighborIndex, R_DimSymbol);
  
  const int minorityRows = INTEGER(minorityDims)[0];
  const int colCount = INTEGER(minorityDims)[1];
  const int neighborRows = INTEGER(neighborDims)[0];
  const int neighborCount = INTEGER(neighborDims)[1];
  
  /* Comprova simetria topologica assegurando que vizinhos batem com matriz de entrada */
  if(neighborRows != minorityRows || XLENGTH(syntheticPerRow) != minorityRows){
    error("Dimensoes inconsistentes para geracao ADASYN");
  }
  
  if(neighborCount < 1){
    error("'minorityNeighborIndex' deve conter ao menos uma coluna");
  }
  
  /* Parametriza arrays fixos de leitura vindos do R */
  const double *minority = REAL(minorityMatrix);
  const int *neighbor = INTEGER(minorityNeighborIndex);
  const int *perRow = INTEGER(syntheticPerRow);
  
  /* Agrega distribuicao almejada para descobrir espaco final demandado da heap */
  R_xlen_t totalSynthetic = 0;
  for(int i = 0; i < minorityRows; ++i){
    if ((i & 8191) == 0) {
      R_CheckUserInterrupt();
    }
    if(perRow[i] < 0){
      error("'syntheticPerRow' nao pode conter valores negativos");
    }
    totalSynthetic += (R_xlen_t) perRow[i];
  }
  
  /* Corta execucoes que excedem capacidade numerica do indice nativo integer */
  if(totalSynthetic > INT_MAX){
    error("Numero de linhas sinteticas excede o limite suportado por matrix R");
  }
  
  /* Pre aloca em memoria de saida protegida todo o limite quantitativo aferido */
  SEXP out = PROTECT(allocMatrix(REALSXP, (int) totalSynthetic, colCount));
  double *synthetic = REAL(out);
  
  /* Importa estado de sementes aleatorias do interpretador garantindo reprodutibilidade */
  GetRNGstate();
  
  /* Anula cursor global usado para rastrear insercoes corretas na matriz expansiva */
  int writeRow = 0;
  for(int i = 0; i < minorityRows; ++i){
    
    /* Escreve e congela estado RNG na memoria antes de permitir chamadas POSIX externas */
    if((i & 255) == 0 && i > 0){
      PutRNGstate();
      R_CheckUserInterrupt();
      GetRNGstate();
    }
    
    const int rowCount = perRow[i];
    
    /* Executa multiplicador de sintese determinado pelo vetor de pesos da linha atual */
    for(int r = 0; r < rowCount; ++r){
      
      /* Utiliza rotina R unif_rand para gerar probabilidade plana e sortear vizinho index */
      int sampledNeighborColumn = (int) floor(unif_rand() * (double) neighborCount);
      
      /* Limita falhas numericas pre condicional caso unif_rand resulte em um literal */
      if(sampledNeighborColumn >= neighborCount){
        sampledNeighborColumn = neighborCount - 1;
      }
      
      /* Rebate sistema de contagem R index 1 para array nativo C index 0 */
      const int selectedNeighborRow = neighbor[i + ((R_xlen_t) sampledNeighborColumn * minorityRows)] - 1;
      
      /* Aciona abort de seguranca caso indices apontem falhas criticas na topologia submetida */
      if(selectedNeighborRow < 0 || selectedNeighborRow >= minorityRows){
        PutRNGstate();
        error("'minorityNeighborIndex' contem indice fora do intervalo");
      }
      
      /* Dispara calculo estocastico de tracao de interpolacao linear delta */
      const double weight = unif_rand();
      
      /* Desenvolve coordenadas pontuais no plano cartesiano feature por feature */
      for(int j = 0; j < colCount; ++j){
        
        /* Resolve ponteiros de salto buscando celulas homologas originais na matriz column-major */
        const double baseValue = minority[i + ((R_xlen_t) j * minorityRows)];
        const double neighborValue = minority[selectedNeighborRow + ((R_xlen_t) j * minorityRows)];
        
        /* Computa diferenca de vetores e impulsiona amostra pelo peso fracionado alocado */
        synthetic[writeRow + ((R_xlen_t) j * totalSynthetic)] = baseValue + weight * (neighborValue - baseValue);
      }
      /* Computa ciclo do gerador de sintese avancando registro da linha destino */
      ++writeRow;
    }
  }
  
  /* Despeja e consolida informacao de saltos RNG na engine R apos fechar iteracao */
  PutRNGstate();
  UNPROTECT(1);
  
  return out;
}


/**
 * @brief Seleciona indices majoritarios NearMiss com menores medias de distancia.
 *
 * Mantem um heap maximo de tamanho retainedMajorityCount para evitar ordenar
 * todos os exemplos majoritarios quando apenas uma fracao pequena sera retida.
 * Empates sao resolvidos pelo menor indice original para espelhar a ordem
 * deterministica usual de R em bases pequenas.
 */
SEXP OU_SelectNearMissMajorityC(SEXP nnDist, SEXP majorityIndex, SEXP retainedMajorityCount){
  require_real_matrix(nnDist, "nnDist");

  if(!isInteger(majorityIndex)){
    error("'majorityIndex' deve ser integer");
  }
  if(!isInteger(retainedMajorityCount) || LENGTH(retainedMajorityCount) != 1){
    error("'retainedMajorityCount' deve ser integer escalar");
  }

  SEXP dims = getAttrib(nnDist, R_DimSymbol);
  const int n = INTEGER(dims)[0];
  const int k = INTEGER(dims)[1];
  const int retain = INTEGER(retainedMajorityCount)[0];

  if(LENGTH(majorityIndex) != n){
    error("'majorityIndex' deve ter comprimento igual a nrow(nnDist)");
  }
  if(retain < 1 || retain > n){
    error("'retainedMajorityCount' deve estar entre 1 e nrow(nnDist)");
  }
  if(k < 1){
    error("'nnDist' deve conter ao menos uma coluna");
  }

  SEXP heapMeanS = PROTECT(allocVector(REALSXP, retain));
  SEXP heapIndexS = PROTECT(allocVector(INTSXP, retain));
  double *heapMean = REAL(heapMeanS);
  int *heapIndex = INTEGER(heapIndexS);
  const double *dist = REAL(nnDist);
  const int *majority = INTEGER(majorityIndex);

  int heapSize = 0;

  #define WORSE(meanA, indexA, meanB, indexB) ((meanA) > (meanB) || ((meanA) == (meanB) && (indexA) > (indexB)))
  #define BETTER(meanA, indexA, meanB, indexB) ((meanA) < (meanB) || ((meanA) == (meanB) && (indexA) < (indexB)))

  for(int i = 0; i < n; ++i){
    if((i & 8191) == 0){
      R_CheckUserInterrupt();
    }

    long double sum = 0.0L;
    for(int j = 0; j < k; ++j){
      const double value = dist[i + ((R_xlen_t) j * n)];
      if(!R_FINITE(value)){
        error("'nnDist' contem distancia ausente ou infinita");
      }
      sum += (long double) value;
    }
    const double mean = (double) (sum / (long double) k);
    const int originalIndex = majority[i];

    if(heapSize < retain){
      int pos = heapSize++;
      heapMean[pos] = mean;
      heapIndex[pos] = originalIndex;
      while(pos > 0){
        int parent = (pos - 1) / 2;
        if(!WORSE(heapMean[pos], heapIndex[pos], heapMean[parent], heapIndex[parent])){
          break;
        }
        double tm = heapMean[parent]; int ti = heapIndex[parent];
        heapMean[parent] = heapMean[pos]; heapIndex[parent] = heapIndex[pos];
        heapMean[pos] = tm; heapIndex[pos] = ti;
        pos = parent;
      }
    }else if(BETTER(mean, originalIndex, heapMean[0], heapIndex[0])){
      heapMean[0] = mean;
      heapIndex[0] = originalIndex;
      int pos = 0;
      while(1){
        int left = 2 * pos + 1;
        int right = left + 1;
        int worst = pos;
        if(left < heapSize && WORSE(heapMean[left], heapIndex[left], heapMean[worst], heapIndex[worst])){
          worst = left;
        }
        if(right < heapSize && WORSE(heapMean[right], heapIndex[right], heapMean[worst], heapIndex[worst])){
          worst = right;
        }
        if(worst == pos){
          break;
        }
        double tm = heapMean[worst]; int ti = heapIndex[worst];
        heapMean[worst] = heapMean[pos]; heapIndex[worst] = heapIndex[pos];
        heapMean[pos] = tm; heapIndex[pos] = ti;
        pos = worst;
      }
    }
  }

  SEXP out = PROTECT(allocVector(INTSXP, retain));
  for(int i = 0; i < retain; ++i){
    INTEGER(out)[i] = heapIndex[i];
  }

  #undef WORSE
  #undef BETTER

  UNPROTECT(3);
  return out;
}

/**
 * @brief Gera sinteticos ADASYN com layout column-friendly e RNG plano.
 *
 * Variante de OU_GenerateSyntheticAdasynC que pre-computa todos os indices
 * de vizinho e pesos de interpolacao em vetores contiguos, e em seguida
 * percorre as colunas no laco externo. Em matrizes column-major isso
 * substitui escritas com stride totalSynthetic por escritas contiguas,
 * melhorando locality de cache em dados de alta dimensionalidade.
 *
 * Resultado numericamente equivalente a OU_GenerateSyntheticAdasynC para
 * o mesmo estado de RNG no momento da chamada.
 */
SEXP OU_GenerateSyntheticAdasynColC(SEXP minorityMatrix, SEXP minorityNeighborIndex, SEXP syntheticPerRow){
  require_real_matrix(minorityMatrix, "minorityMatrix");

  if(!isInteger(minorityNeighborIndex) || !isMatrix(minorityNeighborIndex)){
    error("'minorityNeighborIndex' deve ser uma matrix integer");
  }
  if(!isInteger(syntheticPerRow)){
    error("'syntheticPerRow' deve ser um vetor integer");
  }

  SEXP minorityDims = getAttrib(minorityMatrix, R_DimSymbol);
  SEXP neighborDims = getAttrib(minorityNeighborIndex, R_DimSymbol);
  const R_xlen_t minorityRows = (R_xlen_t) INTEGER(minorityDims)[0];
  const R_xlen_t colCount = (R_xlen_t) INTEGER(minorityDims)[1];
  const R_xlen_t neighborRows = (R_xlen_t) INTEGER(neighborDims)[0];
  const R_xlen_t neighborCount = (R_xlen_t) INTEGER(neighborDims)[1];

  if(neighborRows != minorityRows || XLENGTH(syntheticPerRow) != minorityRows){
    error("Dimensoes inconsistentes para geracao ADASYN");
  }
  if(neighborCount < 1){
    error("'minorityNeighborIndex' deve conter ao menos uma coluna");
  }

  const double *minority = REAL(minorityMatrix);
  const int *neighbor = INTEGER(minorityNeighborIndex);
  const int *perRow = INTEGER(syntheticPerRow);

  /* Calcula total de sinteticos para alocar a matriz final em um unico bloco. */
  R_xlen_t totalSynthetic = 0;
  for(R_xlen_t i = 0; i < minorityRows; ++i){
    if(perRow[i] < 0){
      error("'syntheticPerRow' nao pode conter valores negativos");
    }
    totalSynthetic += (R_xlen_t) perRow[i];
  }
  if(totalSynthetic > INT_MAX){
    error("Numero de linhas sinteticas excede o limite suportado por matrix R");
  }
  if(colCount > INT_MAX){
    error("Numero de colunas excede o limite suportado por matrix R");
  }

  SEXP out = PROTECT(allocMatrix(REALSXP, (int) totalSynthetic, (int) colCount));
  double *synthetic = REAL(out);

  if(totalSynthetic == 0){
    UNPROTECT(1);
    return out;
  }

  /* Pre-aloca vetores temporarios em scratch nao-protegido (R_alloc): para
   * cada sintetico s armazenamos a linha base, a linha vizinha e o peso de
   * interpolacao. Pre-resolver isso permite o laco externo iterar por
   * coluna sem precisar do estado de RNG. */
  int *baseRow = (int *) R_alloc((size_t) totalSynthetic, sizeof(int));
  int *nbrRow  = (int *) R_alloc((size_t) totalSynthetic, sizeof(int));
  double *wgt  = (double *) R_alloc((size_t) totalSynthetic, sizeof(double));

  GetRNGstate();
  R_xlen_t s = 0;
  for(R_xlen_t i = 0; i < minorityRows; ++i){
    if((i & 8191) == 0){
      R_CheckUserInterrupt();
    }
    const int rowCount = perRow[i];
    for(int r = 0; r < rowCount; ++r){
      int sampledCol = (int) floor(unif_rand() * (double) neighborCount);
      if(sampledCol >= (int) neighborCount){
        sampledCol = (int) neighborCount - 1;
      }
      const int selectedNeighborRow = neighbor[i + (R_xlen_t) sampledCol * minorityRows] - 1;
      if(selectedNeighborRow < 0 || selectedNeighborRow >= (int) minorityRows){
        PutRNGstate();
        error("'minorityNeighborIndex' contem indice fora do intervalo");
      }
      baseRow[s] = (int) i;
      nbrRow[s] = selectedNeighborRow;
      wgt[s] = unif_rand();
      ++s;
    }
  }
  PutRNGstate();

  /* Laco externo por coluna: para cada coluna j, percorre todos os sinteticos
   * escrevendo contiguamente na matriz de saida (column-major). Leituras de
   * minority[base + j*minorityRows] tambem sao contiguas para j fixo. */
  for(R_xlen_t j = 0; j < colCount; ++j){
    if((j & 15) == 0){
      R_CheckUserInterrupt();
    }
    const R_xlen_t colOffsetMin = j * minorityRows;
    const R_xlen_t colOffsetSyn = j * totalSynthetic;
    for(R_xlen_t t = 0; t < totalSynthetic; ++t){
      const double bv = minority[baseRow[t] + colOffsetMin];
      const double nv = minority[nbrRow[t] + colOffsetMin];
      synthetic[t + colOffsetSyn] = bv + wgt[t] * (nv - bv);
    }
  }

  UNPROTECT(1);
  return out;
}

/**
 * @brief Remove o proprio ponto de uma matriz de indices KNN em C.
 *
 * Recebe a matriz nbr (n x k_plus) de indices retornada por uma consulta KNN
 * em que query == data, o vetor selfIndex de n indices (a propria linha de
 * cada query) e o numero desejado k de vizinhos validos. Devolve uma matriz
 * n x k onde a coluna que contem self foi removida. Quando self nao aparece
 * em uma linha, mantem os primeiros k vizinhos retornados. Erra com mensagem
 * explicita se sobrar menos do que k candidatos validos.
 */
SEXP OU_DropSelfNeighborC(SEXP nbr, SEXP selfIndex, SEXP desiredK){
  if(!isInteger(nbr) || !isMatrix(nbr)){
    error("'nbr' deve ser uma matrix integer");
  }
  if(!isInteger(selfIndex)){
    error("'selfIndex' deve ser integer");
  }
  if(!isInteger(desiredK) || LENGTH(desiredK) != 1){
    error("'desiredK' deve ser integer escalar");
  }

  SEXP dims = getAttrib(nbr, R_DimSymbol);
  const R_xlen_t n = (R_xlen_t) INTEGER(dims)[0];
  const R_xlen_t kPlus = (R_xlen_t) INTEGER(dims)[1];
  const int k = INTEGER(desiredK)[0];

  if(XLENGTH(selfIndex) != n){
    error("'selfIndex' deve ter comprimento igual a nrow(nbr)");
  }
  if(k < 1){
    error("'desiredK' deve ser >= 1");
  }
  if((R_xlen_t) k > kPlus){
    error("Nao foi possivel remover o proprio ponto mantendo vizinhos suficientes");
  }
  if(n > INT_MAX){
    error("'nbr' excede o limite suportado por allocMatrix");
  }

  SEXP out = PROTECT(allocMatrix(INTSXP, (int) n, k));
  int *outp = INTEGER(out);
  const int *src = INTEGER(nbr);

  for(R_xlen_t i = 0; i < n; ++i){
    if((i & 8191) == 0){
      R_CheckUserInterrupt();
    }
    const int self = INTEGER(selfIndex)[i];
    R_xlen_t written = 0;
    R_xlen_t j = 0;
    /* Percorre as kPlus colunas em ordem; copia ate completar k, pulando
     * a coluna que contem self e NA_INTEGER. */
    while(written < k && j < kPlus){
      const int v = src[i + j * n];
      if(v != NA_INTEGER && v != self){
        outp[i + written * n] = v;
        ++written;
      }
      ++j;
    }
    if(written < k){
      UNPROTECT(1);
      error("Nao foi possivel remover o proprio ponto mantendo vizinhos suficientes");
    }
  }

  UNPROTECT(1);
  return out;
}

/**
 * @brief Top-k brute force KNN para distancia euclidiana usando BLAS.
 *
 * O prototipo de dgemm vem de R_ext/BLAS.h.
 *
 * Calcula D^2[i,j] = ||q_i||^2 + ||r_j||^2 - 2 q_i . r_j em blocos de
 * consulta, via dgemm para o termo cruzado, e usa um max-heap de tamanho
 * k para extrair os k vizinhos com menor distancia para cada query.
 *
 * Entradas:
 *   data:  matriz double n_ref x p (referencia)
 *   query: matriz double n_query x p (consultas)
 *   k:     escalar integer (numero de vizinhos solicitados)
 *
 * Saidas: list(nn.index = integer matrix n_query x k (1-based),
 *              nn.dist  = double matrix n_query x k (distancia euclidiana))
 */
SEXP OU_BruteForceKnnC(SEXP data, SEXP query, SEXP kSEXP){
  require_real_matrix(data, "data");
  require_real_matrix(query, "query");
  if(!isInteger(kSEXP) || LENGTH(kSEXP) != 1){
    error("'k' deve ser integer escalar");
  }

  SEXP dDims = getAttrib(data, R_DimSymbol);
  SEXP qDims = getAttrib(query, R_DimSymbol);
  const R_xlen_t nRef = (R_xlen_t) INTEGER(dDims)[0];
  const R_xlen_t p1 = (R_xlen_t) INTEGER(dDims)[1];
  const R_xlen_t nQuery = (R_xlen_t) INTEGER(qDims)[0];
  const R_xlen_t p2 = (R_xlen_t) INTEGER(qDims)[1];
  const int k = INTEGER(kSEXP)[0];

  if(p1 != p2){
    error("'data' e 'query' devem ter o mesmo numero de colunas");
  }
  if(k < 1){
    error("'k' deve ser >= 1");
  }
  if((R_xlen_t) k > nRef){
    error("'k' nao pode exceder o numero de linhas de 'data'");
  }
  if(nQuery > INT_MAX){
    error("'query' excede o limite suportado por allocMatrix");
  }

  const double *Rref = REAL(data);
  const double *Qref = REAL(query);

  /* Normas ao quadrado por linha, pre-calculadas. */
  double *qNorm2 = (double *) R_alloc((size_t) nQuery, sizeof(double));
  double *rNorm2 = (double *) R_alloc((size_t) nRef, sizeof(double));

  for(R_xlen_t i = 0; i < nQuery; ++i){
    long double s = 0.0L;
    for(R_xlen_t c = 0; c < p1; ++c){
      const double v = Qref[i + c * nQuery];
      s += (long double) v * (long double) v;
    }
    qNorm2[i] = (double) s;
  }
  for(R_xlen_t i = 0; i < nRef; ++i){
    long double s = 0.0L;
    for(R_xlen_t c = 0; c < p1; ++c){
      const double v = Rref[i + c * nRef];
      s += (long double) v * (long double) v;
    }
    rNorm2[i] = (double) s;
  }

  /* Aloca buffers de saida. */
  SEXP outIdx = PROTECT(allocMatrix(INTSXP, (int) nQuery, k));
  SEXP outDst = PROTECT(allocMatrix(REALSXP, (int) nQuery, k));
  int *idxOut = INTEGER(outIdx);
  double *dstOut = REAL(outDst);

  /* Estrategia em blocos para limitar a memoria do produto cruzado a ~64MB. */
  const R_xlen_t maxCells = 8 * 1024 * 1024; /* 8M doubles = 64MB */
  R_xlen_t blockQ = maxCells / (nRef > 0 ? nRef : 1);
  if(blockQ < 1) blockQ = 1;
  if(blockQ > nQuery) blockQ = nQuery;
  if(blockQ > INT_MAX) blockQ = INT_MAX;

  double *cross = (double *) R_alloc((size_t) (blockQ * nRef), sizeof(double));
  /* Buffers para heap: par (dist, idx). */
  double *heapD = (double *) R_alloc((size_t) k, sizeof(double));
  int *heapI = (int *) R_alloc((size_t) k, sizeof(int));

  /* Itera blocos de queries. */
  for(R_xlen_t qStart = 0; qStart < nQuery; qStart += blockQ){
    R_xlen_t qEnd = qStart + blockQ;
    if(qEnd > nQuery) qEnd = nQuery;
    const R_xlen_t curB = qEnd - qStart;
    R_CheckUserInterrupt();

    /* C = Q[qStart:qEnd, ] %*% t(R) : curB x nRef, em column-major. */
    const double alpha = 1.0;
    const double beta = 0.0;
    const int mInt = (int) curB;
    const int nInt = (int) nRef;
    const int kkInt = (int) p1;
    const int ldaInt = (int) nQuery;       /* leading dim de Q (toda matrix) */
    const int ldbInt = (int) nRef;         /* leading dim de R (toda matrix) */
    const int ldcInt = (int) curB;
    /* dgemm calcula C[i,j] = sum_c Q[qStart+i, c] * R[j, c]. Como ambas as
     * matrizes estao em column-major e queremos Q_block @ t(R), usamos
     * transa='N', transb='T' com Q_block apontando para Qref + qStart. */
    F77_CALL(dgemm)(
      "N", "T",
      &mInt, &nInt, &kkInt,
      &alpha,
      Qref + qStart, &ldaInt,
      Rref, &ldbInt,
      &beta,
      cross, &ldcInt FCONE FCONE
    );

    /* Para cada query do bloco, monta heap-top-k de menor distancia. */
    for(R_xlen_t bi = 0; bi < curB; ++bi){
      const R_xlen_t qIdx = qStart + bi;
      const double qn2 = qNorm2[qIdx];
      int heapSize = 0;

      for(R_xlen_t j = 0; j < nRef; ++j){
        /* D^2 = ||q||^2 + ||r||^2 - 2 q.r */
        double d2 = qn2 + rNorm2[j] - 2.0 * cross[bi + j * curB];
        if(d2 < 0.0) d2 = 0.0;  /* clamp para evitar sqrt de negativo por arredondamento */

        if(heapSize < k){
          /* Push: insere no fim e sobe (max-heap por d2). */
          int pos = heapSize++;
          heapD[pos] = d2;
          heapI[pos] = (int) (j + 1);  /* 1-based */
          while(pos > 0){
            int par = (pos - 1) / 2;
            if(heapD[par] >= heapD[pos]) break;
            double td = heapD[par]; int ti = heapI[par];
            heapD[par] = heapD[pos]; heapI[par] = heapI[pos];
            heapD[pos] = td; heapI[pos] = ti;
            pos = par;
          }
        }else if(d2 < heapD[0]){
          /* Substitui raiz e desce. */
          heapD[0] = d2;
          heapI[0] = (int) (j + 1);
          int pos = 0;
          while(1){
            int left = 2 * pos + 1, right = left + 1, worst = pos;
            if(left < heapSize && heapD[left] > heapD[worst]) worst = left;
            if(right < heapSize && heapD[right] > heapD[worst]) worst = right;
            if(worst == pos) break;
            double td = heapD[worst]; int ti = heapI[worst];
            heapD[worst] = heapD[pos]; heapI[worst] = heapI[pos];
            heapD[pos] = td; heapI[pos] = ti;
            pos = worst;
          }
        }
      }

      /* Ordena o heap por distancia crescente (heap sort por seletor). */
      /* Estrategia: extrai-min repetidamente trocando com o final. */
      /* Como o heap atual e max-heap, podemos ordenar in-place produzindo
       * um array em ordem crescente. */
      int sz = heapSize;
      while(sz > 1){
        /* Move raiz (maior) para o final. */
        double td = heapD[0]; int ti = heapI[0];
        heapD[0] = heapD[sz - 1]; heapI[0] = heapI[sz - 1];
        heapD[sz - 1] = td; heapI[sz - 1] = ti;
        --sz;
        /* Restaura heap no prefixo de tamanho sz. */
        int pos = 0;
        while(1){
          int left = 2 * pos + 1, right = left + 1, worst = pos;
          if(left < sz && heapD[left] > heapD[worst]) worst = left;
          if(right < sz && heapD[right] > heapD[worst]) worst = right;
          if(worst == pos) break;
          double td2 = heapD[worst]; int ti2 = heapI[worst];
          heapD[worst] = heapD[pos]; heapI[worst] = heapI[pos];
          heapD[pos] = td2; heapI[pos] = ti2;
          pos = worst;
        }
      }

      /* Escreve resultados na linha qIdx das matrizes de saida (column-major,
       * leading dim = nQuery). Converte distancia^2 para distancia. */
      for(int kk = 0; kk < heapSize; ++kk){
        idxOut[qIdx + (R_xlen_t) kk * nQuery] = heapI[kk];
        dstOut[qIdx + (R_xlen_t) kk * nQuery] = sqrt(heapD[kk]);
      }
    }
  }

  /* Empacota resultado como list(nn.index=..., nn.dist=...) compativel com FNN. */
  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("nn.index"));
  SET_STRING_ELT(names, 1, mkChar("nn.dist"));
  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, outIdx);
  SET_VECTOR_ELT(out, 1, outDst);
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(4);
  return out;
}

/**
 * @brief Mapeia as rotinas internas em C para as chamadas .Call dinamicas provenientes do R.
 */
static const R_CallMethodDef CallEntries[] = {
  {"OU_CheckUserInterruptC", (DL_FUNC) &OU_CheckUserInterruptC, 0},
  {"OU_ComputeZScoreParamsC", (DL_FUNC) &OU_ComputeZScoreParamsC, 1},
  {"OU_ApplyZScoreC", (DL_FUNC) &OU_ApplyZScoreC, 4},
  {"OU_GenerateSyntheticAdasynC", (DL_FUNC) &OU_GenerateSyntheticAdasynC, 3},
  {"OU_GenerateSyntheticAdasynColC", (DL_FUNC) &OU_GenerateSyntheticAdasynColC, 3},
  {"OU_SelectNearMissMajorityC", (DL_FUNC) &OU_SelectNearMissMajorityC, 3},
  {"OU_DropSelfNeighborC", (DL_FUNC) &OU_DropSelfNeighborC, 3},
  {"OU_BruteForceKnnC", (DL_FUNC) &OU_BruteForceKnnC, 3},
  {NULL, NULL, 0}
};

/**
 * @brief Rotina de inicializacao requerida pelo R ao carregar a biblioteca dinamica.
 *
 * Registra o modulo 'instenginer' no interpretador e bloqueia a analise dinamica
 * de simbolos externos para mitigar falhas de resolucao em tempo de execucao.
 *
 * @param dll Referencia referencial de ligacao gerada pela maquina virtual do R.
 */
void attribute_visible R_init_instenginer(DllInfo *dll){
  /* Registra lista de acesso conectando strings R com ponteiros reais no binario */
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  /* Tranca verificacoes dinâmicas reduzindo sobrecarga na busca do kernel */
  R_useDynamicSymbols(dll, FALSE);
  /* Encerra configuracoes assegurando resolucoes declaradas explicitamente */
  R_forceSymbols(dll, TRUE);
}

