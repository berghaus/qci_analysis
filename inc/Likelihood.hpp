#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

#include <map>
#include <vector>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>

#include "Prediction.hpp"
#include "TestStatMonitor.hpp"
#include "PseudoExperiment.hpp"

class TH1;
class TFitterMinuit;
class Prediction;
class TestStatMonitor;

//! Implementation of binned Poisson likelihood
/*! \class Neg2LogLikelihood_FCN
 *
 * binned Poisson likelihood comparing the Experiment given to the constructor to the prediction (also from the
 * constructor) at the given value of the prediction parameter to.  Used in the analysis through the likelihood ratio
 * defined below.
 *
 */
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

//! Interface for likelihood ratio
/*! \class Neg2LogLikelihoodRatio
 *
 * Interface for likelihood ratio implementation.
 *
 */
class Neg2LogLikelihoodRatio {
public:
  Neg2LogLikelihoodRatio();
  Neg2LogLikelihoodRatio( const Experiment* data, const Prediction* pdf, const double& alpha = 0. );
  virtual ~Neg2LogLikelihoodRatio();

  virtual double operator()( const std::vector< double >& ) = 0;
  virtual double Up() const;

  virtual void data( const Experiment* );
  virtual void pdf( const Prediction* );

  virtual const Experiment* data() const;
  virtual const Prediction* pdf() const;
  virtual Neg2LogLikelihood_FCN numerator() const;
  virtual Neg2LogLikelihood_FCN denominator() const;
  virtual void accept( TestStatMonitor& );

protected:
  const Experiment* _data;
  const Prediction* _pdf;

  Neg2LogLikelihood_FCN _numerator;
  Neg2LogLikelihood_FCN _denominator;

};

//! Implementation of maximum likelihood ratio as described in the PDG
/*! \class Neg2LogMaximumLikelihoodRatio
 *
 * The denominator of this likelihood ratio is maximised using MINUIT.
 *
 */
class Neg2LogMaximumLikelihoodRatio: public Neg2LogLikelihoodRatio {

public:

  Neg2LogMaximumLikelihoodRatio();

  // assume ownership of TH1 but not the PDF
  Neg2LogMaximumLikelihoodRatio( const Experiment* data, const Prediction* pdf, const double& alpha = 0. );
  virtual ~Neg2LogMaximumLikelihoodRatio();

  virtual double operator()( const std::vector< double >& );

  virtual void data( const Experiment* );
  virtual void pdf( const Prediction* );

private:

  void init();
  TestStatMonitor * _inCaseShitMonitor; // monitor likelihood denominator fails fit

};


//! Implementation of simple likelihood ratio used by ATLAS and CMS before
/*! \class Neg2LogSimpleLikelihoodRatio
 *
 * The denominator of this likelihood ratio is evaluated against the QCD prediction.
 *
 */
class Neg2LogSimpleLikelihoodRatio: public Neg2LogLikelihoodRatio {

public:

  Neg2LogSimpleLikelihoodRatio();

  // assume ownership of TH1 but not the PDF
  Neg2LogSimpleLikelihoodRatio( const Experiment* data, const Prediction* pdf, const double& alpha = 0. );
  virtual ~Neg2LogSimpleLikelihoodRatio();

  virtual double operator()( const std::vector< double >& );

};


#endif // LIKELIHOOD_HPP
