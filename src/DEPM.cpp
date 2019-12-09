#define TMB_LIB_INIT R_init_DEPM
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(y);
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER(Sigma0);
  PARAMETER(Sigma1);

  int n = x.size();
  vector<Type> F_pred(n);
  vector<Type> Sig_pred(n);

  for(int i = 0; i < n; i++){
    F_pred(i) =  alpha * pow(x(i),beta) ;
    Sig_pred(i) = Sigma0 * pow(F_pred(i),Sigma1) ;
  }

  // Type MLL = 0;
  //
  // for(int i = 0; i < n; i++){
  //   MLL += (log(Sig_pred(i)) + pow(F_pred(i) - y(i),2))/(2*pow(Sig_pred(i),2))    ;
  // }
  //

  Type f;
  f = -sum(dnorm(y,F_pred,Sig_pred, true));

  ADREPORT(F_pred);


  return f;
}
