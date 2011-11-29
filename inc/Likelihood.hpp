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
class PDF;

class Likelihood_FCN: public ROOT::Minuit2::FCNBase {

public:
  Likelihood_FCN();

  // does not assume ownership of TH1 or PDF
  Likelihood_FCN( const Experiment*, const PDF*, const double alpha = 0. );
  virtual ~Likelihood_FCN();

  double operator()() const;
  double operator()( const std::vector< double >& ) const;

  double Up() const;

  void data( const Experiment* );
  void pdf( const PDF* );

  const Experiment* data() const;
  const PDF* pdf() const;
  bool isMinimized() const;

  double Minimize();
  double Minimize( ROOT::Minuit2::MnUserParameters& );

  std::vector< double > pars();

  void accept( TestStatMonitor& );

private:

  const Experiment* _data;
  const PDF* _pdf;

  bool _isMinimized;
  ROOT::Minuit2::MnUserParameters _pars;
  ROOT::Minuit2::MnUserParameterState _parsState;

};

class LikelihoodRatio {

public:

  LikelihoodRatio();

  // assume ownership of TH1 but not the PDF
  LikelihoodRatio( const Experiment* data, const PDF* pdf, const double& alpha = 0. );
  virtual ~LikelihoodRatio();

  double operator()( const std::vector< double >& );

  double Up() const;

  void data( const Experiment* );
  void pdf( const PDF* );

  const Experiment* data() const;
  const PDF* pdf() const;
  Likelihood_FCN numerator() const;   // fit over nuisance parameters
  Likelihood_FCN denominator() const; // global fit

  void accept( TestStatMonitor& );

private:

  void init();

  const Experiment* _data;
  const PDF* _pdf;
  Likelihood_FCN _numerator;   // fit over nuisance parameters
  Likelihood_FCN _denominator; // global fit

};

#endif // LIKELIHOOD_HPP
