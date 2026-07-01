suppressPackageStartupMessages({
  library(sbyadanear)
  library(tibble)
})

ReadIntegerEnv <- function(envName, defaultValue){
  envValue <- Sys.getenv(envName, unset = "")
  
  if(identical(envValue, "")){
    return(defaultValue)
  }
  
  parsedValue <- suppressWarnings(as.integer(envValue))
  
  if(is.na(parsedValue) || parsedValue < 1L){
    return(defaultValue)
  }
  
  return(parsedValue)
}

GetInstanceId <- function(){
  envInstanceId <- Sys.getenv("INSTANCE_ID", unset = "")
  
  if(!identical(envInstanceId, "")){
    return(envInstanceId)
  }
  
  return(paste0(Sys.info()[["nodename"]], "_", Sys.getpid()))
}

CreateTestData <- function(){
  set.seed(20260515L)
  
  majorityCount <- 180L
  minorityCount <- 20L
  predictorCount <- 8L
  
  majorityMatrix <- matrix(
    data = rnorm(majorityCount * predictorCount, mean = 0, sd = 1),
    nrow = majorityCount,
    ncol = predictorCount
  )
  
  minorityMatrix <- matrix(
    data = rnorm(minorityCount * predictorCount, mean = 1.8, sd = 0.9),
    nrow = minorityCount,
    ncol = predictorCount
  )
  
  majorityMatrix[, 7L] <- rbinom(majorityCount, size = 1L, prob = 0.25)
  minorityMatrix[, 7L] <- rbinom(minorityCount, size = 1L, prob = 0.65)
  
  majorityMatrix[, 8L] <- rpois(majorityCount, lambda = 1)
  minorityMatrix[, 8L] <- rpois(minorityCount, lambda = 3)
  
  predictorMatrix <- rbind(majorityMatrix, minorityMatrix)
  
  targetVector <- c(
    rep("majority", majorityCount),
    rep("minority", minorityCount)
  )
  
  rowOrder <- sample.int(length(targetVector))
  
  testData <- data.frame(
    predictorMatrix[rowOrder, , drop = FALSE],
    TARGET = factor(targetVector[rowOrder]),
    stringsAsFactors = FALSE
  )
  
  names(testData) <- c(
    sprintf("x%02d", seq_len(predictorCount)),
    "TARGET"
  )
  
  return(as_tibble(testData))
}

ExtractBalancedData <- function(resultObject){
  if(is.data.frame(resultObject)){
    return(as_tibble(resultObject))
  }
  
  candidateNames <- c(
    "sby_balanced_data",
    "sby_final_data",
    "sby_data",
    "sby_under_data",
    "sby_over_data"
  )
  
  for(candidateName in candidateNames){
    if(!is.null(resultObject[[candidateName]]) && is.data.frame(resultObject[[candidateName]])){
      return(as_tibble(resultObject[[candidateName]]))
    }
  }
  
  stop("Nao foi possivel localizar o tibble balanceado no objeto retornado")
}

ValidateBalancedData <- function(balancedData){
  if(!is.data.frame(balancedData)){
    stop("O retorno balanceado nao e data.frame ou tibble")
  }
  
  if(!"TARGET" %in% names(balancedData)){
    stop("A coluna TARGET nao foi encontrada no retorno")
  }
  
  targetTable <- table(balancedData$TARGET)
  
  if(length(targetTable) != 2L){
    stop("O retorno nao contem exatamente duas classes em TARGET")
  }
  
  if(any(targetTable < 1L)){
    stop("O retorno contem classe sem observacoes")
  }
  
  return(targetTable)
}

RunCase <- function(methodName, testData, workerCount, repetitionIndex){
  seedValue <- 100000L + repetitionIndex
  
  if(identical(methodName, "adasyn")){
    resultObject <- sby_adasyn(
      sby_formula = TARGET ~ .,
      sby_data = testData,
      sby_ratio_over = 0.50,
      sby_knn_over_k = 5L,
      sby_seed = seedValue,
      sby_audit = FALSE,
      sby_knn_engine = "auto",
      sby_knn_algorithm = "auto",
      sby_knn_distance_metric = "euclidean",
      sby_knn_workers = workerCount
    )
  }
  
  if(identical(methodName, "nearmiss")){
    resultObject <- sby_nearmiss(
      sby_formula = TARGET ~ .,
      sby_data = testData,
      sby_ratio_under = 0.80,
      sby_knn_under_k = 5L,
      sby_seed = seedValue,
      sby_audit = FALSE,
      sby_knn_engine = "auto",
      sby_knn_algorithm = "auto",
      sby_knn_distance_metric = "euclidean",
      sby_knn_workers = workerCount
    )
  }
  
  if(identical(methodName, "adanear")){
    resultObject <- sby_adanear(
      sby_formula = TARGET ~ .,
      sby_data = testData,
      sby_ratio_over = 0.50,
      sby_ratio_under = 0.80,
      sby_knn_over_k = 5L,
      sby_knn_under_k = 5L,
      sby_seed = seedValue,
      sby_audit = FALSE,
      sby_knn_engine = "auto",
      sby_knn_algorithm = "auto",
      sby_knn_distance_metric = "euclidean",
      sby_knn_workers = workerCount
    )
  }
  
  balancedData <- ExtractBalancedData(resultObject)
  targetTable <- ValidateBalancedData(balancedData)
  
  return(data.frame(
    methodName = methodName,
    repetitionIndex = repetitionIndex,
    rowsOutput = nrow(balancedData),
    majorityCount = as.integer(targetTable[["majority"]]),
    minorityCount = as.integer(targetTable[["minority"]]),
    status = "ok",
    stringsAsFactors = FALSE
  ))
}



CollectMklDiagnostics <- function(){
  cfg <- sbyadanear:::sby_resolve_oneapi_mkl()
  data.frame(
    methodName = "mkl_diagnostics",
    repetitionIndex = 0L,
    rowsOutput = NA_integer_,
    majorityCount = NA_integer_,
    minorityCount = NA_integer_,
    status = if(isTRUE(cfg$enabled)) "ok" else "ok",
    elapsedSeconds = 0,
    message = paste0(
      "enabled=", as.character(cfg$enabled),
      "; threads=", as.character(cfg$threads),
      "; OMP_NUM_THREADS=", Sys.getenv("OMP_NUM_THREADS", unset = ""),
      "; MKL_NUM_THREADS=", Sys.getenv("MKL_NUM_THREADS", unset = "")
    ),
    stringsAsFactors = FALSE
  )
}

Main <- function(){
  instanceId <- GetInstanceId()
  repetitionCount <- ReadIntegerEnv("TEST_REPETITIONS", 2L)
  workerCount <- ReadIntegerEnv("TEST_KNN_WORKERS", 2L)
  outputDir <- Sys.getenv("TEST_OUTPUT_DIR", unset = "test-results")
  
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  
  testData <- CreateTestData()
  methodNames <- c("adasyn", "nearmiss", "adanear")
  
  resultList <- list()
  resultIndex <- 1L
  
  for(repetitionIndex in seq_len(repetitionCount)){
    for(methodName in methodNames){
      startTime <- proc.time()[["elapsed"]]
      
      runResult <- tryCatch(
        expr = {
          caseResult <- RunCase(
            methodName = methodName,
            testData = testData,
            workerCount = workerCount,
            repetitionIndex = repetitionIndex
          )
          
          caseResult$elapsedSeconds <- round(proc.time()[["elapsed"]] - startTime, 6)
          caseResult$message <- ""
          caseResult
        },
        error = function(errorCondition){
          data.frame(
            methodName = methodName,
            repetitionIndex = repetitionIndex,
            rowsOutput = NA_integer_,
            majorityCount = NA_integer_,
            minorityCount = NA_integer_,
            status = "error",
            elapsedSeconds = round(proc.time()[["elapsed"]] - startTime, 6),
            message = conditionMessage(errorCondition),
            stringsAsFactors = FALSE
          )
        }
      )
      
      resultList[[resultIndex]] <- runResult
      resultIndex <- resultIndex + 1L
    }
  }
  
  resultList[[resultIndex]] <- CollectMklDiagnostics()
  resultData <- do.call(rbind, resultList)
  resultData$instanceId <- instanceId
  resultData$workerCount <- workerCount
  
  outputPath <- file.path(
    outputDir,
    paste0("sbyadanear_tests_", gsub("[^A-Za-z0-9_]+", "_", instanceId), ".csv")
  )
  
  write.csv(
    x = resultData,
    file = outputPath,
    row.names = FALSE
  )
  
  print(resultData)
  
  failureCount <- sum(resultData$status != "ok")
  
  if(failureCount > 0L){
    quit(status = 1L, save = "no")
  }
  
  quit(status = 0L, save = "no")
}

Main()
