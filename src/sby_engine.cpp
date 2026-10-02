#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <mkl.h>
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>
#include <dlfcn.h>

#ifndef _OPENMP
#error "sbyadanear requires Intel oneAPI OpenMP; compile with icpx -qopenmp."
#endif
#ifndef __INTEL_LLVM_COMPILER
#error "sbyadanear requires the Intel oneAPI icpx compiler."
#endif

namespace {
using Clock = std::chrono::steady_clock;
struct ThreadScope {
  int omp_before, dynamic_before, mkl_before;
  explicit ThreadScope(int n): omp_before(omp_get_max_threads()),
      dynamic_before(omp_get_dynamic()), mkl_before(mkl_set_num_threads_local(n)) {
    omp_set_dynamic(0);
    omp_set_num_threads(n);
  }
  ~ThreadScope() {
    mkl_set_num_threads_local(mkl_before);
    omp_set_num_threads(omp_before);
    omp_set_dynamic(dynamic_before);
  }
};
double rss() {
  std::ifstream f("/proc/self/status");
  std::string key, rest;
  while(f >> key) {
    if(key == "VmRSS:") { double kb; f >> kb; return kb * 1024.; }
    std::getline(f, rest);
  }
  return NA_REAL;
}
double peak_rss() {
  struct rusage r;
  return getrusage(RUSAGE_SELF, &r) == 0 ? r.ru_maxrss * 1024. : NA_REAL;
}
struct Meter {
  std::vector<std::string> stages;
  std::vector<double> seconds, memory_before, memory_after, peak;
  std::vector<int> observed;
  const bool enabled;
  explicit Meter(bool active): enabled(active) {}
  Clock::time_point t;
  double before;
  void start() { if(enabled) { t = Clock::now(); before = rss(); } }
  void finish(const std::string& name, int threads) {
    if(!enabled) return;
    stages.push_back(name);
    seconds.push_back(std::chrono::duration<double>(Clock::now()-t).count());
    memory_before.push_back(before); memory_after.push_back(rss());
    peak.push_back(peak_rss()); observed.push_back(threads);
  }
  Rcpp::DataFrame frame(int threads) const {
    const int n = stages.size();
    return Rcpp::DataFrame::create(
      Rcpp::_["stage"] = stages, Rcpp::_["elapsed_seconds"] = seconds,
      Rcpp::_["rss_before_bytes"] = memory_before,
      Rcpp::_["rss_after_bytes"] = memory_after,
      Rcpp::_["process_peak_rss_bytes"] = peak,
      Rcpp::_["omp_observed_threads"] = observed,
      Rcpp::_["thread_limit"] = Rcpp::IntegerVector(n, threads),
      Rcpp::_["mkl_configured_threads"] = Rcpp::IntegerVector(n, mkl_get_max_threads()),
      Rcpp::_["mkl_observed_threads"] = Rcpp::IntegerVector(n, NA_INTEGER));
  }
};
struct Matrix {
  int n, p;
  std::vector<double> a;
  // Matrices become immutable before neighbor queries. Cache row norms only
  // within this call; the same extended-precision accumulation is reused.
  mutable std::vector<double> norms;
  Matrix(int nr, int nc): n(nr), p(nc), a(static_cast<size_t>(nr)*nc) {}
  double& at(int i,int j) { return a[static_cast<size_t>(i)*p+j]; }
  double at(int i,int j) const { return a[static_cast<size_t>(i)*p+j]; }
};
Matrix subset(const Matrix& x, const std::vector<int>& ids, int& observed) {
  Matrix out(ids.size(), x.p);
  if(!x.norms.empty()) out.norms.resize(ids.size());
#pragma omp parallel for schedule(static) reduction(max:observed)
  for(int i=0;i<out.n;++i) {
    observed = std::max(observed, omp_get_num_threads());
    if(!out.norms.empty()) out.norms[i]=x.norms[ids[i]];
    std::copy_n(&x.a[static_cast<size_t>(ids[i])*x.p], x.p,
                &out.a[static_cast<size_t>(i)*x.p]);
  }
  return out;
}
void ensure_norms(const Matrix& x, int& observed) {
  if(!x.norms.empty()) return;
  x.norms.resize(x.n);
#pragma omp parallel for schedule(static) reduction(max:observed)
  for(int i=0;i<x.n;++i) {
    observed=std::max(observed,omp_get_num_threads());
    long double sum=0;
    for(int j=0;j<x.p;++j) sum+=static_cast<long double>(x.at(i,j))*x.at(i,j);
    x.norms[i]=static_cast<double>(sum);
  }
}
struct Neighbor { double d; int index; };
// DGEMM is outside OpenMP regions: each BLAS call can use the same Intel
// thread budget without nested teams. Only a bounded query tile is allocated.
std::vector<Neighbor> knn(const Matrix& ref, const Matrix& query, int k,
                         bool farthest, const std::vector<int>& self,
                         int& observed, int tile) {
  std::vector<Neighbor> out(static_cast<size_t>(query.n)*k);
  ensure_norms(ref,observed);
  ensure_norms(query,observed);
  const auto& norm_ref=ref.norms;
  const auto& norm_query=query.norms;
  // Allocate on the calling thread: allocation failures must unwind through
  // Rcpp and ThreadScope, rather than escape an OpenMP worker.
  std::vector<std::vector<Neighbor>> scratch(omp_get_max_threads());
  for(auto& row:scratch) row.reserve(k);
  // DGEMM beta=0 overwrites the active tile; reuse storage between blocks.
  std::vector<double> dot(static_cast<size_t>(std::min(tile,query.n))*ref.n);
  for(int first=0;first<query.n;first+=tile) {
    Rcpp::checkUserInterrupt();
    int count=std::min(tile,query.n-first);
    int invalid_distance=0;
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,count,ref.n,ref.p,
      1.,&query.a[static_cast<size_t>(first)*query.p],query.p,
      ref.a.data(),ref.p,0.,dot.data(),ref.n);
#pragma omp parallel for schedule(static) reduction(max:observed,invalid_distance)
    for(int r=0;r<count;++r) {
      observed = std::max(observed, omp_get_num_threads());
      int qi=first+r;
      auto& row=scratch[omp_get_thread_num()];
      row.clear();
      auto less=[farthest](const Neighbor& a,const Neighbor& b) {
        return a.d==b.d ? a.index<b.index : (farthest ? a.d>b.d : a.d<b.d);
      };
      // A bounded heap retains only K exact distances. Its root is the worst
      // retained neighbor, including the original-index tie-breaker.
      for(int i=0;i<ref.n;++i) {
        if(!self.empty() && i==self[qi]) continue;
        double cross=dot[static_cast<size_t>(r)*ref.n+i];
        double d2=norm_query[qi]+norm_ref[i]-2.*cross;
        double bound=64.*std::numeric_limits<double>::epsilon()*
          (norm_query[qi]+norm_ref[i]+2.*std::abs(cross));
        // Recompute near cancellation; coincident rows remain legitimate
        // neighbors, and only the exact self index is excluded.
        double distance;
        if(!std::isfinite(d2) || d2<=bound) {
          long double sum=0;
          for(int j=0;j<ref.p;++j) {
            long double d=static_cast<long double>(query.at(qi,j))-ref.at(i,j);
            sum += d*d;
          }
          d2=static_cast<double>(sum);
          // Square roots in extended precision preserve finite distances when
          // the squared value overflows or underflows double storage.
          distance=(!std::isfinite(d2) || (d2==0. && sum>0.)) ?
            static_cast<double>(std::sqrt(sum)) : std::sqrt(std::max(0.,d2));
        } else distance=std::sqrt(std::max(0.,d2));
        if(!std::isfinite(distance)) invalid_distance=1;
        const Neighbor candidate{distance,i};
        if(static_cast<int>(row.size())<k) {
          row.push_back(candidate);
          if(static_cast<int>(row.size())==k) std::make_heap(row.begin(),row.end(),less);
        } else if(less(candidate,row.front())) {
          std::pop_heap(row.begin(),row.end(),less);
          row.back()=candidate;
          std::push_heap(row.begin(),row.end(),less);
        }
      }
      std::sort_heap(row.begin(),row.end(),less);
      std::copy_n(row.begin(),k,out.begin()+static_cast<size_t>(qi)*k);
    }
    if(invalid_distance) Rcpp::stop("Euclidean distances exceed double precision; revise the supplied scaling information.");
  }
  return out;
}
Rcpp::IntegerVector index1(const std::vector<int>& ids) {
  Rcpp::IntegerVector out(ids.size());
  for(size_t i=0;i<ids.size();++i) out[i]=ids[i]+1;
  return out;
}
Rcpp::NumericMatrix rmatrix(const Matrix& x) {
  Rcpp::NumericMatrix out(x.n,x.p);
  for(int i=0;i<x.n;++i) for(int j=0;j<x.p;++j) out(i,j)=x.at(i,j);
  return out;
}
Rcpp::IntegerVector affinity() {
  cpu_set_t set; CPU_ZERO(&set);
  std::vector<int> cpus;
  if(sched_getaffinity(0,sizeof(set),&set)==0)
    for(int i=0;i<CPU_SETSIZE;++i) if(CPU_ISSET(i,&set)) cpus.push_back(i);
  return Rcpp::wrap(cpus);
}
std::string library_path(void *symbol) {
  Dl_info info;
  return dladdr(symbol,&info) && info.dli_fname ? info.dli_fname : "unknown";
}
}

extern "C" SEXP sby_pipeline_cpp(SEXP xs, SEXP ys, SEXP os) {
  BEGIN_RCPP
  Rcpp::NumericMatrix input(xs);
  Rcpp::IntegerVector y(ys);
  Rcpp::List opt(os);
  int n=input.nrow(),p=input.ncol();
  int threads=Rcpp::as<int>(opt["threads"]), rare=Rcpp::as<int>(opt["rare"]);
  int g=Rcpp::as<int>(opt["synthetics"]), model=Rcpp::as<int>(opt["model"]);
  int ko=Rcpp::as<int>(opt["ko"]),ku=Rcpp::as<int>(opt["ku"]);
  int km=Rcpp::as<int>(opt["km"]),tile=Rcpp::as<int>(opt["tile"]);
  bool detailed=Rcpp::as<bool>(opt["audit"]);
  bool uniform=Rcpp::as<bool>(opt["uniform"]);
  double under=Rcpp::as<double>(opt["under"]);
  if(n<2 || p<1 || y.size()!=n || threads<1 || g<0 ||
     model<1 || model>3 || ko<1 || ku<1 || km<1 || tile<1 ||
     !std::isfinite(under) || under<0)
    Rcpp::stop("Invalid internal pipeline parameters.");
  ThreadScope scope(threads);
  Meter meter(detailed);
  meter.start();
  Matrix raw(n,p),z(n,p);
  std::vector<int> minority,majority;
  for(int i=0;i<n;++i) {
    if(y[i]==NA_INTEGER) Rcpp::stop("Missing class label.");
    (y[i]==rare ? minority : majority).push_back(i);
    for(int j=0;j<p;++j) {
      if(!std::isfinite(input(i,j))) Rcpp::stop("Predictors must be finite.");
      raw.at(i,j)=input(i,j);
    }
  }
  if(minority.empty() || majority.empty()) Rcpp::stop("Exactly two observed classes are required.");
  if(static_cast<long long>(n)+g>std::numeric_limits<int>::max())
    Rcpp::stop("Output row count exceeds the supported integer range.");
  int observed=0;
  Rcpp::NumericVector center(p),scale(p);
  // Optional externally supplied scale is validated in R. Otherwise it is
  // computed ONCE on the original data, with population denominator n.
  bool external=opt.containsElementNamed("centers");
  if(external) {
    center=Rcpp::clone(Rcpp::NumericVector(opt["centers"]));
    scale=Rcpp::clone(Rcpp::NumericVector(opt["scales"]));
  }
  std::vector<int> invalid(p,0);
#pragma omp parallel for schedule(static) reduction(max:observed)
  for(int j=0;j<p;++j) {
    observed=std::max(observed,omp_get_num_threads());
    long double mean=0,variance=0;
    if(!external) {
      for(int i=0;i<n;++i) mean+=raw.at(i,j);
      mean/=n;
      for(int i=0;i<n;++i) {
        long double d=static_cast<long double>(raw.at(i,j))-mean;
        variance+=d*d;
      }
      center[j]=static_cast<double>(mean);
      scale[j]=static_cast<double>(std::sqrt(variance/n));
      if(scale[j]==0) scale[j]=1.;
    }
    if(!std::isfinite(center[j]) || !std::isfinite(scale[j]) || scale[j]<=0) {
      invalid[j]=1; continue;
    }
    for(int i=0;i<n;++i) {
      z.at(i,j)=static_cast<double>((static_cast<long double>(raw.at(i,j))-center[j])/scale[j]);
      if(!std::isfinite(z.at(i,j))) invalid[j]=1;
    }
  }
  if(std::accumulate(invalid.begin(),invalid.end(),0))
    Rcpp::stop("Initial scaling is not representable in double precision.");
  meter.finish("initial_scaling",observed);
  meter.start(); observed=0;
  const int nm=minority.size();
  Matrix mz=subset(z,minority,observed);
  std::vector<int> hits(nm,0),quota(nm,0);
  bool used_uniform=false;
  int effective_ko=std::min(ko,n-1), effective_min_k=std::min(ko,nm-1);
  if(g>0) {
    if(nm<2) Rcpp::stop("ADASYN requer ao menos dois registros raros originais.");
    auto mixed=knn(z,mz,effective_ko,false,minority,observed,tile);
    long long total=0;
    for(int i=0;i<nm;++i) {
      for(int j=0;j<effective_ko;++j)
        if(y[mixed[static_cast<size_t>(i)*effective_ko+j].index]!=rare) ++hits[i];
      total+=hits[i];
    }
    if(total==0 && !uniform)
      Rcpp::stop("ADASYN: todas as dificuldades sao zero. Use sby_adasyn_zero_difficulty = 'uniform' para solicitar a extensao documentada.");
    used_uniform=total==0;
    std::vector<long double> fraction(nm);
    long long allocated=0;
    for(int i=0;i<nm;++i) {
      long double exact=static_cast<long double>(g)*(used_uniform ? 1.L/nm : static_cast<long double>(hits[i])/total);
      quota[i]=static_cast<int>(std::floor(exact));
      fraction[i]=exact-quota[i]; allocated+=quota[i];
    }
    std::vector<int> order(nm); std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b) {
      return fraction[a]==fraction[b] ? minority[a]<minority[b] : fraction[a]>fraction[b];
    });
    long long remaining=static_cast<long long>(g)-allocated;
    if(remaining<0 || remaining>nm) Rcpp::stop("Invalid ADASYN integer allocation.");
    for(int i=0;i<remaining;++i) ++quota[order[i]];
  }
  meter.finish("adasyn_difficulty_and_quotas",observed);
  meter.start(); observed=0;
  Matrix syn(g,p),sz(g,p);
  std::vector<int> parent(g),partner(g);
  std::vector<double> lambda(g);
  if(g>0) {
    std::vector<int> self(nm); std::iota(self.begin(),self.end(),0);
    auto nearby=knn(mz,mz,effective_min_k,false,self,observed,tile);
    // RNG calls remain on the R thread; every synthetic uses one lambda for
    // all columns. RNGkind is fixed/restored by the scoped R caller.
    Rcpp::RNGScope rng;
    int s=0;
    for(int i=0;i<nm;++i) {
      if(i%1024==0) Rcpp::checkUserInterrupt();
      for(int j=0;j<quota[i];++j) {
      int choice=std::min(effective_min_k-1,static_cast<int>(R::runif(0.,1.)*effective_min_k));
      parent[s]=minority[i];
      partner[s]=minority[nearby[static_cast<size_t>(i)*effective_min_k+choice].index];
      lambda[s]=R::runif(0.,1.); ++s;
      }
    }
#pragma omp parallel for schedule(static) reduction(max:observed)
    for(int i=0;i<g;++i) {
      observed=std::max(observed,omp_get_num_threads());
      for(int j=0;j<p;++j) {
        sz.at(i,j)=std::lerp(z.at(parent[i],j),z.at(partner[i],j),lambda[i]);
        // Affine inverse expressed with original parents avoids cancellation
        // when rare values are tiny compared with the initial center.
        syn.at(i,j)=std::lerp(raw.at(parent[i],j),raw.at(partner[i],j),lambda[i]);
      }
    }
  }
  meter.finish("adasyn_interpolation_and_inverse",observed);
  meter.start(); observed=0;
  std::vector<int> selected=majority,candidates=majority;
  std::vector<double> scores(majority.size(),NA_REAL);
  int target=majority.size(),effective_ku=std::min(ku,nm+g),effective_km=std::min(km,static_cast<int>(majority.size()));
  bool shortage=false;
  // Cap as double BEFORE conversion: even a finite ratio near DBL_MAX is safe.
  if(under>0) target=static_cast<int>(std::min(static_cast<double>(majority.size()),
                                  std::floor((static_cast<double>(nm)+g)*under)));
  if(under>0 && target<1) Rcpp::stop("NearMiss: retencao arredonda para zero; aumente sby_nearmiss_ratio ou use zero para desativar.");
  const bool nearmiss_executed=under>0 && target<static_cast<int>(majority.size());
  if(nearmiss_executed) {
    Matrix rare_z(nm+g,p);
    std::copy(mz.a.begin(),mz.a.end(),rare_z.a.begin());
    std::copy(sz.a.begin(),sz.a.end(),rare_z.a.begin()+mz.a.size());
    Matrix major_z=subset(z,majority,observed);
    std::vector<bool> eligible(majority.size(),model!=3);
    if(model==3) {
      auto pre=knn(major_z,rare_z,effective_km,false,{},observed,tile);
      for(const auto& v:pre) eligible[v.index]=true;
    }
    candidates.clear();
    std::vector<int> local;
    for(size_t i=0;i<majority.size();++i) if(eligible[i]) {
      candidates.push_back(majority[i]); local.push_back(i);
    }
    Matrix q=subset(major_z,local,observed);
    auto distances=knn(rare_z,q,effective_ku,model==2,{},observed,tile);
#pragma omp parallel for schedule(static) reduction(max:observed)
    for(int i=0;i<static_cast<int>(local.size());++i) {
      observed=std::max(observed,omp_get_num_threads());
      long double sum=0;
      for(int j=0;j<effective_ku;++j) sum+=distances[i*effective_ku+j].d;
      scores[local[i]]=static_cast<double>(sum/effective_ku);
    }
    shortage=static_cast<int>(local.size())<target;
    target=std::min(target,static_cast<int>(local.size()));
    // The tie-breaker makes the order total; sorting the retained prefix
    // selects exactly the same rows without ordering discarded candidates.
    std::partial_sort(local.begin(),local.begin()+target,local.end(),[&](int a,int b) {
      return scores[a]==scores[b] ? majority[a]<majority[b] :
             (model==3 ? scores[a]>scores[b] : scores[a]<scores[b]);
    });
    selected.resize(target);
    for(int i=0;i<target;++i) selected[i]=majority[local[i]];
    std::sort(selected.begin(),selected.end());
  }
  meter.finish("nearmiss_selection",observed);
  Rcpp::List result=Rcpp::List::create(
    Rcpp::_["synthetic"] = rmatrix(syn),
    Rcpp::_["minority_indices"] = index1(minority),
    Rcpp::_["majority_indices"] = index1(selected),
    Rcpp::_["centers"] = center, Rcpp::_["scales"] = scale,
    Rcpp::_["uniform_fallback"] = used_uniform,
    Rcpp::_["adasyn_executed"] = g>0,
    Rcpp::_["nearmiss_executed"] = nearmiss_executed,
    Rcpp::_["candidate_shortage"] = shortage,
    Rcpp::_["effective_k"] = Rcpp::IntegerVector::create(effective_ko,effective_min_k,effective_ku,effective_km),
    Rcpp::_["telemetry"] = meter.frame(threads));
  if(detailed) {
    result["parents"]=index1(parent); result["partners"]=index1(partner);
    result["lambda"]=Rcpp::wrap(lambda); result["synthetic_z"]=rmatrix(sz);
    result["difficulty_hits"]= g>0 ? Rcpp::wrap(hits) : Rcpp::wrap(Rcpp::IntegerVector(nm,NA_INTEGER));
    result["quotas"]=Rcpp::wrap(quota);
    result["candidates"]=index1(candidates); result["scores"]=Rcpp::wrap(scores);
  }
  return result;
  END_RCPP
}
extern "C" SEXP sby_runtime_cpp() {
  BEGIN_RCPP
  char version[256]; mkl_get_version_string(version,sizeof(version));
  std::string omp_path=library_path(reinterpret_cast<void*>(&omp_get_max_threads));
  if(omp_path.find("libiomp5")==std::string::npos)
    Rcpp::stop("Intel libiomp5 was not resolved; check for conflicting OpenMP runtimes.");
  return Rcpp::List::create(Rcpp::_["openmp"]=true,
    Rcpp::_["mkl_linked"]=true,Rcpp::_["compiler"]="Intel oneAPI icpx",
    Rcpp::_["mkl_version"]=version,Rcpp::_["openmp_library"]=omp_path,
    Rcpp::_["mkl_library"]=library_path(reinterpret_cast<void*>(&mkl_get_version_string)),
    Rcpp::_["affinity_cpus"]=affinity(),
    Rcpp::_["openmp_max_threads"]=omp_get_max_threads(),
    Rcpp::_["openmp_thread_limit"]=omp_get_thread_limit(),
    Rcpp::_["openmp_dynamic"]=omp_get_dynamic(),
    Rcpp::_["mkl_max_threads"]=mkl_get_max_threads());
  END_RCPP
}
extern "C" SEXP sby_thread_probe_cpp(SEXP ts,SEXP ns) {
  BEGIN_RCPP
  int t=Rcpp::as<int>(ts),n=Rcpp::as<int>(ns),observed=0;
  if(t<1 || n<2 || n>8192) Rcpp::stop("Probe: threads >= 1 and 2 <= n <= 8192 required.");
  ThreadScope scope(t);
#pragma omp parallel reduction(max:observed)
  { observed=std::max(observed,omp_get_num_threads()); }
  std::vector<double> a(static_cast<size_t>(n)*n,1.),b(a.size(),2.),c(a.size());
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.,a.data(),n,b.data(),n,0.,c.data(),n);
  return Rcpp::List::create(Rcpp::_["omp_observed_threads"]=observed,
      Rcpp::_["mkl_configured_threads"]=mkl_get_max_threads(),
      Rcpp::_["checksum"]=c[0]);
  END_RCPP
}
static const R_CallMethodDef calls[] = {
  {"sby_pipeline_cpp",(DL_FUNC)&sby_pipeline_cpp,3},
  {"sby_runtime_cpp",(DL_FUNC)&sby_runtime_cpp,0},
  {"sby_thread_probe_cpp",(DL_FUNC)&sby_thread_probe_cpp,2},
  {NULL,NULL,0}
};
extern "C" void attribute_visible R_init_sbyadanear(DllInfo *dll) {
  R_registerRoutines(dll,NULL,calls,NULL,NULL);
  R_useDynamicSymbols(dll,static_cast<Rboolean>(0));
  R_forceSymbols(dll,static_cast<Rboolean>(0));
}
