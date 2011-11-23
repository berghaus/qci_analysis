#include "PValueTest.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/checked_delete.hpp>

#include <TCanvas.h>
#include <TFitterMinuit.h>
#include <TH1.h>

#include "PDF.hpp"
#include "Likelihood.hpp"

#define foreach BOOST_FOREACH
using namespace std;
using boost::format;

PValueTest::PValueTest( const double alpha, const vector< LikelihoodRatio* >& lambdas ) :
    _alpha( alpha ),
    _lambdas( lambdas ) {
  init();
}

void PValueTest::init() {

  string hName = str( format( "Likelihood_FCN-alpha%2.1e" ) % _alpha );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(), "", 1000, 0., -1. );
  string xTitle = str( format( "-2ln#lambda(%2.0f TeV )" ) % _alpha );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );

  vector< double > par( 1, _alpha );
  _testStats.reserve( _lambdas.size() );
  foreach( LikelihoodRatio* lambda, _lambdas ) {
    LikelihoodRatio& l = *lambda;
    _testStats.push_back( l( par ) );
  }
  sort( _testStats.begin(), _testStats.end() );

}

PValueTest::PValueTest() :
    _alpha( 0 ) {
}

PValueTest::~PValueTest() {
}

double PValueTest::operator()( LikelihoodRatio& lambda ) {

  vector< double > par( 1, _alpha );
  double dataLL = lambda( par );
  vector< double >::iterator itr = lower_bound( _testStats.begin(), _testStats.end(), dataLL );
  int nOver = distance( itr, _testStats.end() );

  return double( nOver ) / double( _testStats.size() );

}

void PValueTest::finalize() {
  foreach( const double& x, _testStats ) {
    _minus2LnLikelihoodDistribution->Fill( x );
  }

  TCanvas* pvc = new TCanvas( "PValueCanvas", "", 500, 500 );
  pvc->cd();
  _minus2LnLikelihoodDistribution->Draw();
  string cName = _minus2LnLikelihoodDistribution->GetName();
  pvc->Print( ( cName + ".png" ).c_str() );
}

