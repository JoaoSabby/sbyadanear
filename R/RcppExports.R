# Native symbol bindings for registered .Call interfaces.
# Symbols are resolved lazily because package installation byte-compiles R code
# before every registered native routine is available through the loaded DLL.

sby_native_symbol_cache <- new.env(parent = emptyenv())

sby_load_native_symbol <- function(sby_symbol_name){
  if(exists(sby_symbol_name, envir = sby_native_symbol_cache, inherits = FALSE)){
    return(get(
      x = sby_symbol_name,
      envir = sby_native_symbol_cache,
      inherits = FALSE
    ))
  }

  sby_symbol <- getNativeSymbolInfo(
    name = sby_symbol_name,
    PACKAGE = "sbyadanear"
  )

  assign(
    x = sby_symbol_name,
    value = sby_symbol,
    envir = sby_native_symbol_cache
  )

  return(sby_symbol)
}

sby_native_symbol_available <- function(sby_symbol_name){
  tryCatch(
    {
      sby_load_native_symbol(sby_symbol_name)
      TRUE
    },
    error = function(sby_error){
      FALSE
    }
  )
}

sby_call_native <- function(sby_symbol_name, ...){
  .Call(
    sby_load_native_symbol(sby_symbol_name),
    ...
  )
}

# Observacao sobre resolucao de simbolos nativos:
# O pacote usa useDynLib(sbyadanear, .registration = TRUE), portanto o R ja
# cria automaticamente, no namespace do pacote, um objeto NativeSymbolInfo para
# cada rotina .Call registrada em R_init_sbyadanear. A resolucao em tempo de
# execucao e feita de forma preguicosa por sby_load_native_symbol("nome"), que
# consulta getNativeSymbolInfo() e mantem um cache em sby_native_symbol_cache.
# Nao criamos aqui binds adicionais por delayedAssign(): alem de redundantes,
# eles colidiam com os objetos gerados pela auto-registracao, gerando avisos
# 'failed to assign RegisteredNativeSymbol ... already defined in namespace' e a
# falha fatal 'no such symbol' durante o lazy loading.
