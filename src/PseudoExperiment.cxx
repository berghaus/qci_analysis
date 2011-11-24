#include "PseudoExperiment.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/checked_delete.hpp>

#include <TCanvas.h>
#include <TH2.h>

using namespace std;
using namespace boost;

PseudoExperimentFactory::PseudoExperimentFactory( const PDF* pdf, const PseudoExperiment* graft, unsigned int seed ) :
    _pdf( pdf ),
    _graft( graft ),
    _random( seed ) {
}

PseudoExperimentFactory::PseudoExperimentFactory( const PDF* pdf, const TH1* graft, unsigned int seed ) :
    _pdf( pdf ),
    _graft( static_cast< const PseudoExperiment* >( graft ) ),
    _random( seed ) {
}

PseudoExperimentFactory::~PseudoExperimentFactory() {

  // can't delete PEs since root does not like deleting TH1's
//  typedef std::map< double, std::vector< PseudoExperiment* > > pes_t;
//  foreach( pes_t::value_type& x, _generated )
//    for_each( x.second.begin(), x.second.end(), checked_deleter< PseudoExperiment >() );

}

PseudoExperiment *
PseudoExperimentFactory::build( const double& alpha ) {

  if ( _nGenerated.find( alpha ) == _nGenerated.end() ) _nGenerated[alpha] = 1;
  else ++_nGenerated[alpha];

  string peName = str( format( "PE_alpha%1.0e_n%1.0f" ) % alpha % _nGenerated[alpha] );
  PseudoExperiment* result = (PseudoExperiment*) _graft->Clone( peName.c_str() );
  result->Reset( "ICES" );

  for( int bin = 1; bin <= _graft->GetNbinsX(); ++bin ) {
    double chi = _graft->GetBinCenter( bin );
    double expectedN = ( *_pdf )( chi, alpha );
    int n = _random.Poisson( expectedN );
    result->SetBinContent( bin, n );
  }

  _generated[alpha].push_back( result );
  return result;
}

vector< PseudoExperiment* > PseudoExperimentFactory::build( const double& alpha, const int& n ) {

  vector< PseudoExperiment* > result( n );

  if ( _generated.find( alpha ) == _generated.end() || _generated.find( alpha )->second.size() < n ) {

    do {
      build( alpha );
    } while ( _generated.find( alpha )->second.size() < n );
    copy( _generated.find( alpha )->second.begin(), _generated.find( alpha )->second.end(), result.begin() );

  } else {

    copy( _generated.find( alpha )->second.begin(), _generated.find( alpha )->second.begin() + n, result.begin() );

  }
  return result;
}

PseudoExperimentFactory::PseudoExperimentFactory() :
    _pdf( 0 ),
    _graft( 0 ) {
}
