#ifndef PSEUDO_EXPERIMENT_HPP
#define PSEUDO_EXPERIMENT_HPP
#include <map>
#include <vector>

#include <TH1.h>
#include <TRandom.h>

#include "PDF.hpp"


typedef TH1D PseudoExperiment;

class PseudoExperimentFactory {

public:
  PseudoExperimentFactory( const PDF* pdf, const PseudoExperiment* graft, unsigned int seed = 65539 );
  PseudoExperimentFactory( const PDF* pdf, const TH1* graft, unsigned int seed = 65539 );
  virtual ~PseudoExperimentFactory();
  std::vector<PseudoExperiment*> build( const double& alpha, const int& n );

private:

  PseudoExperimentFactory();
  const PDF *              _pdf;
  std::map<double,TH1*>    _profiles;
  const PseudoExperiment * _graft;

  TRandom _random;
  std::map<double,unsigned int>      _nGenerated;
  std::map<double,std::vector<PseudoExperiment*> > _generated;

  PseudoExperiment * build( const double& alpha );

};


#endif // PSEUDO_EXPERIMENT_HPP
