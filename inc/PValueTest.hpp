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
  PValueTest( const double, const LikelihoodRatio_FCN& , std::vector<PseudoExperiment*> );
  virtual ~PValueTest();

  double operator() ( const TH1* );
  void init();

private:
  PValueTest();
  const double                   _alpha; // alpha used for PEs
  LikelihoodRatio_FCN            _testStat;
  std::vector<PseudoExperiment*> _pes;

  TH1 * _minus2LnLikelihoodDistribution;
  std::vector<double> _testStats;

};


#endif // PVALUE_TEST_HPP
