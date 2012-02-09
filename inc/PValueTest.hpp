#ifndef PVALUE_TEST_HPP
#define PVALUE_TEST_HPP
#include <iostream>
#include <vector>
#include "PseudoExperiment.hpp"
#include "Likelihood.hpp"

class TFitterMinuit;
class TH1;
class Prediction;

class PValueTest {

public:
  PValueTest();
  PValueTest( const double, const std::vector< Neg2LogLikelihoodRatio* >& );
  PValueTest( const double&, const int&, PseudoExperimentFactory& );
  virtual ~PValueTest();

  double operator()( Neg2LogLikelihoodRatio& );
  double operator()( const double& );
  void finalize();
  double alpha() const;
  void alpha( const double& );

  friend std::ostream& operator<<( std::ostream&, const PValueTest& );
  friend std::istream& operator>>( std::istream&, PValueTest& );

private:

  void init();
  void init( const int&, PseudoExperimentFactory& ); //<- construct likelihood distribution from PEs
  void clear(); //<- purge all information in PValueTest


  double _alpha; // alpha used for PEs
  std::vector< Neg2LogLikelihoodRatio* > _lambdas;
  double _dataLLR;

  TH1 * _minus2LnLikelihoodDistribution;
  std::vector< double > _testStats;

};

#endif // PVALUE_TEST_HPP
