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
  
  /* Extrai metadados do objeto para capturar quantidades de linhas e colunas */
  SEXP dims = getAttrib(xMatrix, R_DimSymbol);
  const int n = INTEGER(dims)[0]; 
  const int p = INTEGER(dims)[1]; 
  
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
  for(int j = 0; j < p; ++j){
    /* Previne sobrecarga no processador verificando interrupcoes apenas a cada 16 ciclos */
    if((j & 15) == 0){
      R_CheckUserInterrupt();
    }
    
    /* Posiciona o ponteiro de leitura no byte zero da coluna atual */
    const double *col = x + ((R_xlen_t) j * n);
    long double sum = 0.0L;
    
    /* Executa primeira passagem agregando valores para derivar o momento inicial absoluto */
    for(int i = 0; i < n; ++i){
      sum += (long double) col[i];
    }
    
    /* Define o parametro de centralizacao e grava no vetor protegido */
    const long double mean = sum / (long double) n;
    mu[j] = (double) mean;
    
    long double ss = 0.0L;
    
    /* Executa segunda passagem acumulando quadrados residuais estritamente sobre residuos centrados */
    for (int i = 0; i < n; ++i) {
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
  const int n = INTEGER(dims)[0];
  const int p = INTEGER(dims)[1];
  
  /* Assegura conformidade de espaco vetorial de colunas x vetores de distribuicao */
  if(XLENGTH(centers) != p || XLENGTH(scales) != p){
    error("'centers' e 'scales' devem ter comprimento igual ao numero de colunas");
  }
  
  /* Cria espelho da matriz alvo alocando novo bloco de memoria homogeneo */
  SEXP out = PROTECT(allocMatrix(REALSXP, n, p));
  
  /* Estabelece vinculacoes de enderecos para varredura C padrao */
  const double *x = REAL(xMatrix);
  const double *mu = REAL(centers);
  const double *sd = REAL(scales);
  double *y = REAL(out);
  
  /* Traduz sinalizador booleano logico do R em inteiro avaliavel em C */
  const int do_reverse = asLogical(reverse);
  
  for(int j = 0; j < p; ++j){
    /* Restringe consultas de hook de sistema operacional */
    if((j & 15) == 0){
      R_CheckUserInterrupt();
    }
    
    /* Pre calcula salto contiguo de bytes para a base do eixo y */
    const R_xlen_t offset = (R_xlen_t) j * n;
    const double center = mu[j];
    const double scale = sd[j];
    
    if(do_reverse){
      /* Restaura posicao real baseando se em distribuicao escalada pre formatada */
      for(int i = 0; i < n; ++i){
        y[offset + i] = x[offset + i] * scale + center;
      }
    } else {
      /* Normaliza vetor realocando ao redor de zero em desvios unificados unitarios */
      for(int i = 0; i < n; ++i){
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
 * @brief Mapeia as rotinas internas em C para as chamadas .Call dinamicas provenientes do R.
 */
static const R_CallMethodDef CallEntries[] = {
  {"OU_CheckUserInterruptC", (DL_FUNC) &OU_CheckUserInterruptC, 0},
  {"OU_ComputeZScoreParamsC", (DL_FUNC) &OU_ComputeZScoreParamsC, 1},
  {"OU_ApplyZScoreC", (DL_FUNC) &OU_ApplyZScoreC, 4},
  {"OU_GenerateSyntheticAdasynC", (DL_FUNC) &OU_GenerateSyntheticAdasynC, 3},
  {"OU_SelectNearMissMajorityC", (DL_FUNC) &OU_SelectNearMissMajorityC, 3},
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

