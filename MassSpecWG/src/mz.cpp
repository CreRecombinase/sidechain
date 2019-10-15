#include <Rcpp.h>



//[[Rcpp::export]]
Rcpp::NumericVector filter_vec(Rcpp::NumericVector target, Rcpp::NumericVector query,const double delta=1e-04){

  std::vector<double> return_vec;
  return_vec.reserve(query.size());
  auto l_b = target.begin();
  auto target_e = target.end();
  for(auto q : query){
    l_b = std::lower_bound(l_b,target_e,q-q*delta);
    if(l_b == target_e)
      break;

    auto u_b = std::lower_bound(l_b,target_e,q+q*delta);
    std::copy(l_b,u_b,std::back_inserter(return_vec));
    if(u_b==target_e){
      break;
    }
  }
  return(Rcpp::wrap(return_vec));
}



//[[Rcpp::export]]
Rcpp::LogicalVector filter_pred_range(Rcpp::NumericVector targetLow, Rcpp::NumericVector targetHigh, Rcpp::NumericVector query,const double delta=1e-04){

  Rcpp::LogicalVector ret(targetLow.size(),FALSE);
  auto q_b = query.begin();
  auto q_e = query.end();
  auto min_ql = *q_b-delta*(*q_b);
  auto max_ql= *q_b+delta*(*q_b);
  const size_t p= targetLow.size();
  double t_tl;
  double t_th;
  for(int i=0; i<p; i++){
    t_tl=targetLow(i);
    t_th=targetHigh(i);
    while (t_tl > max_ql){
      if (q_b == q_e){
        return ret;
      }
      q_b++;
      min_ql = *q_b-delta*(*q_b);
      max_ql= *q_b+delta*(*q_b);
    }
    if(((t_tl <= max_ql)) && ((min_ql <= t_th))){
      ret(i)= TRUE;
    }
  }
  return(ret);
}



//[[Rcpp::export]]
Rcpp::LogicalVector filter_pred(Rcpp::NumericVector target, Rcpp::NumericVector query,const double delta=1e-04){
  Rcpp::LogicalVector ret(target.size(),FALSE);
  auto q_b = query.begin();
  auto q_e = query.end();
  auto min_ql = *q_b-delta*(*q_b);
  auto max_ql= *q_b+delta*(*q_b);
  const size_t p= target.size();
  double t_t;
  for(int i=0; i<p; i++){
    t_t=target(i);
    while (t_t > max_ql){
      if (q_b == q_e){
        return ret;
      }
      q_b++;
      min_ql = *q_b-delta*(*q_b);
      max_ql= *q_b+delta*(*q_b);
    }
    if(t_t <= max_ql && (t_t > min_ql)){
      ret(i)= TRUE;
    }
  }
  return(ret);
}
