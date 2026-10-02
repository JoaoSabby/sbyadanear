fixture <- function() {
  data.frame(x=c(.2,.7,1.2,1.8, .25,.5,.8,1.1,1.4,1.7,2.2,2.8,3.3,4.1),
    z=c(.4,1.1,.6,1.7, .3,1.2,.7,1.5,.2,1.4,2.1,2.4,3.5,4.4),
    binary=c(0,1,0,1,0,1,0,1,0,1,0,1,0,1),
    integer=as.integer(c(-2,0,2,4,-2,-1,0,1,2,3,4,5,6,7)),
    whole_double=as.double(0:13),
    constant=rep(1L,14),
    target=ordered(c(rep("rare",4),rep("major",10)),levels=c("major","unused","rare")))
}
values <- function(x) lapply(x,identity)
oracle_nearmiss <- function(z, rare, major, k, model, m, retained) {
  d <- outer(seq_len(nrow(major)),seq_len(nrow(rare)),Vectorize(function(i,j)
    sqrt(sum((major[i,]-rare[j,])^2))))
  candidates <- seq_len(nrow(major))
  if (model==3L) candidates <- sort(unique(unlist(lapply(seq_len(nrow(rare)),function(j)
    order(d[,j],seq_len(nrow(major)))[seq_len(min(m,nrow(major)))]))))
  score <- apply(d,1,function(v) mean(sort(v,decreasing=model==2L)[seq_len(min(k,length(v)))]))
  chosen <- candidates[order(if(model==3L) -score[candidates] else score[candidates],candidates)]
  sort(chosen[seq_len(min(retained,length(chosen)))])
}
sampler <- function(d=fixture(),...) sby_adanear_hpc(d,target~.,sby_seed=17L,
  sby_adasyn_k=3L,sby_nearmiss_k=2L,sby_config_max_threads=2L,...)
