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
  PValueTest( const double, const std::vector<LikelihoodRatio*>& );
  virtual ~PValueTest();

  double operator() ( LikelihoodRatio& );
  void init();
  void finalize();

private:
  PValueTest();
  const double                  _alpha; // alpha used for PEs
  std::vector<LikelihoodRatio*> _lambdas;

  TH1 * _minus2LnLikelihoodDistribution;
  std::vector<double> _testStats;

};


#endif // PVALUE_TEST_HPP
