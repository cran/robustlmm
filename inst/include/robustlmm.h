#ifndef __robustlmm_h__
#define __robustlmm_h__

#include <R_ext/Applic.h>
#include <Rcpp.h>

class Integration {
private:
  int neval_, ier_, limit_, lenw_, last_;
  double epsabs_, epsrel_, result_, noBound_, abserr_;
  
  int* iwork_;
  double* work_;
  
public: 
  Integration();
  ~Integration(); 
  
  /* typedef void integr_fn(double *x, int n, void *ex); */
  double integrateNInfInf(integr_fn f, void *ex);
  double integrateAInf(integr_fn f, void *ex, double* bound);
  double integrateNInfB(integr_fn f, void *ex, double* bound);
  double integrateAB(integr_fn f, void *ex, double* a, double* b);
  
  int getNeval();
  int getIer();
  int getLast();
  double getAbserr();
  
private:
  void dqagi(integr_fn f, void *ex, double* bound, int inf);
  void dqags(integr_fn f, void *ex, double* a, double* b);
  void checkIer();
};

class PsiFunction {
public:
  PsiFunction();
  
  virtual const std::string name() const;
  virtual const std::string show() const;
  virtual void chgDefaults(Rcpp::NumericVector tDefs);
  virtual Rcpp::NumericVector tDefs() const;
  virtual const std::string showDefaults() const;
  
  virtual const double rhoFun(const double x);
  virtual const double psiFun(const double x);
  virtual const double wgtFun(const double x);
  virtual const double DpsiFun(const double x);
  virtual const double DwgtFun(const double x);
  const double psi2Fun(const double x);
  
  virtual const double Erho();
  virtual const double Epsi2();
  virtual const double EDpsi();
  
  virtual ~PsiFunction();
};

typedef const double (PsiFunction::*Fptr)(const double);

class PsiFunctionNumIntExp : public PsiFunction {
public:
  PsiFunctionNumIntExp();
  
  const std::string name() const;
  void chgDefaults(Rcpp::NumericVector tDefs);
  
  virtual const double Erho();
  virtual const double Epsi2();
  virtual const double EDpsi();
  
  ~PsiFunctionNumIntExp();
  
private:
  double Erho_;
  double Epsi2_;
  double EDpsi_;
  Integration integration_;
  
  void reset();
  
  const double computeErho();
  const double computeEpsi2();
  const double computeEDpsi();
  double integrate(Fptr fptr);
};

#endif // __robustlmm_h__