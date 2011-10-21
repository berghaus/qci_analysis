#ifndef PVALUE_TEST_HPP
#define PVALUE_TEST_HPP

#include "PseudoExperiment.hpp"


class PValueTest {

public:
  PValueTest();
  virtual ~PValueTest();

  double operator() ( const double& );

private:
  PseudoExperimentFactory* _factory;

};


#endif // PVALUE_TEST_HPP
