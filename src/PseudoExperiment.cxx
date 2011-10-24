#include "PseudoExperiment.hpp"
#include <iostream>
#include <string>

#include <boost/format.hpp>

#include <TCanvas.h>
#include <TH2.h>

using namespace std;

using boost::format;

PseudoExperimentFactory::PseudoExperimentFactory( const PDF* pdf,
						  const PseudoExperiment* graft,
						  unsigned int seed )
 : _pdf       ( pdf )
 , _graft     ( graft )
 , _random    ( seed ) {
}


PseudoExperimentFactory::PseudoExperimentFactory( const PDF* pdf,
						  const TH1* graft,
						  unsigned int seed ) 
 : _pdf       ( pdf )
 , _graft     ( static_cast<const PseudoExperiment*>(graft) )
 , _random    ( seed ) {
}


PseudoExperimentFactory::~PseudoExperimentFactory() {
}


PseudoExperiment *
PseudoExperimentFactory::build( const double& alpha ) {

  if ( _nGenerated.find( alpha ) == _nGenerated.end() ) _nGenerated[alpha] = 1;
  else ++_nGenerated[alpha];

  string peName = str( format("PE_alpha%1.0e_n%1.0f") % alpha % _nGenerated[alpha] );
  PseudoExperiment* result = (PseudoExperiment*)_graft->Clone( peName.c_str() );
  result->Reset("ICES");

  for ( int bin = 1; bin <= _graft->GetNbinsX(); ++bin ) {
    double chi       = _graft->GetBinCenter( bin );    
    double expectedN = (*_pdf)( chi, alpha );
    int n            = _random.Poisson( expectedN );
    result->SetBinContent( bin, n );
  }

  return result;
}


PseudoExperimentFactory::PseudoExperimentFactory()
  : _pdf  (0)
  , _graft(0) {
}
