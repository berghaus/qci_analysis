#ifndef PVALUE_TEST_HPP
#define PVALUE_TEST_HPP
#include <iostream>
#include <string>
#include <vector>
#include "PseudoExperiment.hpp"
#include "Likelihood.hpp"

class TFitterMinuit;
class TH1;
class Prediction;

template< typename Likelihood >
class PValueTest {

public:
  PValueTest();
  PValueTest( const double, const std::vector< Likelihood* >& );

  virtual ~PValueTest();

  double operator()( Likelihood& );
  double operator()( const double& );
  void finalize( const std::string& dir = "./" );
  double alpha() const;
  void alpha( const double& );

  template< typename L >
  friend std::ostream& operator<<( std::ostream&, const PValueTest<L>& );

  template< typename L >
  friend std::istream& operator>>( std::istream&, PValueTest<L>& );

private:

  void init();
  void clear(); //<- purge all information in PValueTest

  double _alpha; // alpha used for PEs
  std::vector< Likelihood* > _lambdas;
  double _dataLLR;

  TH1 * _minus2LnLikelihoodDistribution;
  std::vector< double > _testStats;

};

#endif // PVALUE_TEST_HPP
