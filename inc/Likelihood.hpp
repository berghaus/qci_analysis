#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

#include <map>
#include <vector>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>

#include "PDF.hpp"
#include "TestStatMonitor.hpp"

class TH1;
class TFitterMinuit;
class PDF;

class Likelihood_FCN : public ROOT::Minuit2::FCNBase {

public:
  Likelihood_FCN();

  // does not assume ownership of TH1 or PDF
  Likelihood_FCN( const TH1*, const PDF* );
  virtual ~Likelihood_FCN();

  double operator() ( ) const;
  double operator() ( const std::vector<double>& ) const;

  double Up() const;


  void data( const TH1* );
  void pdf ( const PDF* );

  const TH1* data() const;
  const PDF* pdf () const;
  bool       isMinimized() const;

  double Minimize();
  double Minimize( ROOT::Minuit2::MnUserParameters& );

  std::vector<double> pars();

  void accept( TestStatMonitor& );

private:

  const TH1* _data;
  const PDF* _pdf;


  bool _isMinimized;
  ROOT::Minuit2::MnUserParameters     _pars;
  ROOT::Minuit2::MnUserParameterState _parsState;

};



class LikelihoodRatio_FCN {

public:

  LikelihoodRatio_FCN();

  // assume ownership of TH1 but not the PDF
  LikelihoodRatio_FCN( const TH1*, const PDF* );
  virtual ~LikelihoodRatio_FCN();

  double operator() ( const std::vector<double>& );

  double Up() const;

  void data( const TH1* );
  void pdf ( const PDF* );

  TH1* data() const;
  const PDF* pdf () const;

  void accept( TestStatMonitor& );

private:

  void init();

  TH1* _data;
  const PDF* _pdf;
  Likelihood_FCN _numerator;   // fit over nuisance parameters
  Likelihood_FCN _denominator; // global fit

};


#endif // LIKELIHOOD_HPP
