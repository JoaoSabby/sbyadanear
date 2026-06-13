\name{sby_hpc_cpu_report}
\alias{sby_hpc_cpu_report}
\title{Relatorio de CPU e compilacao do motor HPC}
\usage{
sby_hpc_cpu_report()
}
\description{
Consolida evidencias para confirmar se o servidor e o pacote estao prontos para
usar o caminho Intel Cascade Lake com AVX-512, AVX2, FMA e OpenMP.
}
\details{
A confirmacao e dividida em runtime_flags, compile_report e binary_scan. Para a
meta do cliente, considere OK quando runtime_ok, compile_ok e openmp_ok forem
TRUE. A varredura do binario por zmm/vfmadd via objdump e uma evidencia adicional
quando objdump esta instalado.
}
\value{
Lista com indicadores logicos e evidencias de runtime, compilacao, ambiente
NUMA/MKL e varredura opcional do binario nativo.
}
