#include "PseudoExperiment.hpp"

#include <TH2.h>

using namespace std;


PseudoExperimentFactory::PseudoExperimentFactory( const PDF* pdf, const PseudoExperiment* graft )
 : _pdf  ( pdf )
 , _graft( graft ) {
}


PseudoExperimentFactory::PseudoExperimentFactory( const PDF* pdf, const TH1* graft ) 
 : _pdf  ( pdf )
 , _graft( static_cast<const PseudoExperiment*>(graft) ) {
}


PseudoExperimentFactory::~PseudoExperimentFactory() {
}


PseudoExperiment *
PseudoExperimentFactory::build( const double& alpha  ) {
  PseudoExperiment* result = (PseudoExperiment*)_graft->Clone("PE");
  return result;
}


PseudoExperimentFactory::PseudoExperimentFactory()
  : _pdf  (0)
  , _graft(0)
{}
