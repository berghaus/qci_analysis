#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

#include <map>
#include <vector>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>

#include "PDF.hpp"
#include "TestStatMonitor.hpp"
#include "PseudoExperiment.hpp"

class TH1;
class TFitterMinuit;
class Prediction;
class TestStatMonitor;

class Neg2LogLikelihood_FCN: public ROOT::Minuit2::FCNBase {

public:
  Neg2LogLikelihood_FCN();

  // does not assume ownership of TH1 or PDF
  Neg2LogLikelihood_FCN( const Experiment*, const Prediction*, const double alpha = 0. );
  virtual ~Neg2LogLikelihood_FCN();

  double operator()() const;
  double operator()( const std::vector< double >& ) const;

  double Up() const;

  void data( const Experiment* );
  void pdf( const Prediction* );

  const Experiment* data() const;
  const Prediction* pdf() const;
  bool isMinimized() const;

  double Minimize();
  double Minimize( ROOT::Minuit2::MnUserParameters& );

  std::vector< double > pars() const;
  void pars( const std::vector< double >& );

  void accept( TestStatMonitor& );

private:

  const Experiment* _data;
  const Prediction* _pdf;

  bool _isMinimized;
  ROOT::Minuit2::MnUserParameters _pars;
  ROOT::Minuit2::MnUserParameterState _parsState;

};

class Neg2LogLikelihoodRatio {

public:

  Neg2LogLikelihoodRatio();

  // assume ownership of TH1 but not the PDF
  Neg2LogLikelihoodRatio( const Experiment* data, const Prediction* pdf, const double& alpha = 0. );
  virtual ~Neg2LogLikelihoodRatio();

  double operator()( const std::vector< double >& );

  double Up() const;

  void data( const Experiment* );
  void pdf( const Prediction* );

  const Experiment* data() const;
  const Prediction* pdf() const;
  Neg2LogLikelihood_FCN numerator() const;   // fit over nuisance parameters
  Neg2LogLikelihood_FCN denominator() const; // global fit

  void accept( TestStatMonitor& );

private:

  void init();

  const Experiment* _data;
  const Prediction* _pdf;
  Neg2LogLikelihood_FCN _numerator;   // fit over nuisance parameters
  Neg2LogLikelihood_FCN _denominator; // global fit
  TestStatMonitor * _inCaseShitMonitor; // monitor likelihood denominator fails fit

};

#endif // LIKELIHOOD_HPP
