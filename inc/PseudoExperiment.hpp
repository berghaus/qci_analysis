#ifndef PSEUDO_EXPERIMENT_HPP
#define PSEUDO_EXPERIMENT_HPP
#include <map>

#include <TH1.h>
#include <TRandom.h>

#include "PDF.hpp"


typedef TH1D PseudoExperiment;

class PseudoExperimentFactory {

public:
  PseudoExperimentFactory( const PDF* pdf, const PseudoExperiment* graft, unsigned int seed = 65539 );
  PseudoExperimentFactory( const PDF* pdf, const TH1* graft, unsigned int seed = 65539 );
  virtual ~PseudoExperimentFactory();
  PseudoExperiment * build( const double& alpha = 0. );

private:

  PseudoExperimentFactory();
  const PDF *              _pdf;
  std::map<double,TH1*>    _profiles;
  const PseudoExperiment * _graft;

  TRandom _random;
  std::map<double,unsigned int> _nGenerated;

};


#endif // PSEUDO_EXPERIMENT_HPP
