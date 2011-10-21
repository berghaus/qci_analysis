#ifndef PSEUDO_EXPERIMENT_HPP
#define PSEUDO_EXPERIMENT_HPP
#include <map>
#include <TH1.h>

#include "PDF.hpp"


typedef TH1D PseudoExperiment;

class PseudoExperimentFactory {

  PseudoExperimentFactory( const PDF*, const PseudoExperiment* );
  virtual ~PseudoExperimentFactory();
  PseudoExperiment build( const double& alpha = 0. );

private:

  PseudoExperimentFactory();
  const PDF *        _pdf;
  std::map<double,TH1*>   _profiles;
  const PseudoExperiment * _graft;

};


#endif // PSEUDO_EXPERIMENT_HPP
