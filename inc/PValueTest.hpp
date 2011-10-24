#ifndef PVALUE_TEST_HPP
#define PVALUE_TEST_HPP
#include <vector>
#include "PseudoExperiment.hpp"
#include "Likelihood.hpp"

class TFitterMinuit;
class TH1;
class PDF;

class PValueTest {

public:
  PValueTest( const double, const LikelihoodRatio& , PseudoExperimentFactory* );
  virtual ~PValueTest();

  double operator() ( const TH1* );
  void init( unsigned int nPE = 10000 );

private:
  PValueTest( );
  const double             _alpha; // alpha used for PEs
  LikelihoodRatio          _testStat;
  PseudoExperimentFactory* _factory;

  TH1 * _minus2LnLikelihoodDistribution;

  unsigned int              _nPE;

};


#endif // PVALUE_TEST_HPP
