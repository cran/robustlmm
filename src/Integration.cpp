#include "Integration.h"

#include<sstream>
template <typename T>
std::string to_string(T value)
{
  std::ostringstream os ;
  os << value ;
  return os.str() ;
}

#define DEBUG(STRING)                                          \
// Rcpp::Rcout << STRING << std::endl;

Integration::Integration() : neval_(0), ier_(0), limit_(100), 
  lenw_(4*limit_), last_(0), 
  epsabs_(std::pow(std::numeric_limits<double>::epsilon(), .5)),
  epsrel_(epsabs_), result_(0), noBound_(NA_REAL), abserr_(0),
  iwork_(Calloc(limit_, int)), work_(Calloc(lenw_, double)) {}

Integration::~Integration() {
  Free(iwork_); Free(work_);
}

double Integration::integrateNInfInf(integr_fn f, void *ex) {
  dqagi(f, ex, &noBound_, 2);
  return result_;
}

double Integration::integrateAInf(integr_fn f, void *ex, double* bound) {
  dqagi(f, ex, bound, 1);
  return result_;
}

double Integration::integrateNInfB(integr_fn f, void *ex, double* bound) {
  dqagi(f, ex, bound, -1);
  return result_;
}

double Integration::integrateAB(integr_fn f, void *ex, double* a, double* b) {
  dqags(f, ex, a, b);
  return result_;
}

int Integration::getNeval() {
  return neval_;
}

int Integration::getIer() {
  return ier_;
}

int Integration::getLast() {
  return last_; 
}

double Integration::getAbserr() {
  return abserr_;
}

void Integration::dqagi(integr_fn f, void *ex, double* bound, int inf) {
  DEBUG("Calling Rdqagi with inf = " << inf << " and bound = " << *bound)
  Rdqagi(f, ex, bound, &inf, &epsabs_, &epsrel_, &result_,
         &abserr_, &neval_, &ier_, &limit_, &lenw_, &last_, 
         iwork_, work_);
  DEBUG(" Result: " << result_)
  checkIer();
}
  
void Integration::dqags(integr_fn f, void *ex, double* a, double* b) {
  DEBUG("Calling Rdqags with a = " << *a << " and b = " << *b)
  Rdqags(f, ex, a, b, &epsabs_, &epsrel_, &result_,
         &abserr_, &neval_, &ier_, &limit_, &lenw_, &last_, 
         iwork_, work_);
  DEBUG(" Result: " << result_)
  checkIer();
}

void Integration::checkIer() {
  if(ier_ > 0 && ier_!=5) 
    Rcpp::warning("integration flag " +  to_string(ier_)); 
}

/*
void Integration::reset() {
  neval_ = 0;
  ier_ = 0;
  last_ = 0;
  result_ = 0.;
  abserr_ = 0.;
} 
*/
