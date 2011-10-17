#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

#include <map>
#include <vector>

#include <Minuit2/FCNBase.h>

#include "PDF.hpp"

class TH1;
class TFitterMinuit;
class PDF;

class Likelihood : public ROOT::Minuit2::FCNBase {

public:
  Likelihood();

  // does not assume ownership of TH1 or PDF
  Likelihood( const TH1*, const PDF* );
  virtual ~Likelihood();

  double operator() ( const std::vector<double>& ) const;

  double Up() const;


  void data( const TH1* );
  void pdf ( const PDF* );

  const TH1* data() const;
  const PDF* pdf () const;
  

private:

  const TH1* _data;
  const PDF* _pdf;

};



class LikelihoodRatio {

public:

  LikelihoodRatio();

  // assume ownership of TH1 but not the PDF
  LikelihoodRatio( TFitterMinuit *, const TH1*, const PDF* );
  virtual ~LikelihoodRatio();

  double operator() ( const std::vector<double>& ) const;

  double Up() const;

  void data( const TH1* );
  void pdf ( const PDF* );

  TH1* data() const;
  const PDF* pdf () const;

private:

  void init();

  TH1* _data;
  const PDF* _pdf;
  Likelihood _numerator;   // fit over nuisance parameters
  double     _denominator; // result of global fit
  Likelihood _denominatorL; // global fit
  TFitterMinuit * _fitter; // does not seem to like getting deleted

};


#endif // LIKELIHOOD_HPP
