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
#include <TLine.h>

#include "PDF.hpp"
#include "Likelihood.hpp"

#define foreach BOOST_FOREACH
using namespace std;
using boost::format;

PValueTest::PValueTest( const double alpha, const vector< Neg2LogLikelihoodRatio* >& lambdas ) :
    _alpha( alpha ),
    _lambdas( lambdas ),
    _dataLLR( 0 ) {
  init();
}

void PValueTest::init() {

  vector< double > par( 1, _alpha );
  _testStats.reserve( _lambdas.size() );
  foreach( Neg2LogLikelihoodRatio* lambda, _lambdas )
  {
    Neg2LogLikelihoodRatio& l = *lambda;
    _testStats.push_back( l( par ) );
  }
  sort( _testStats.begin(), _testStats.end() );

}

PValueTest::PValueTest() :
    _alpha( 0 ),
    _dataLLR( 0 ) {
}

PValueTest::~PValueTest() {
}

double PValueTest::operator()( Neg2LogLikelihoodRatio& lambda ) {

  vector< double > par( 1, _alpha );
  _dataLLR = lambda( par );
  vector< double >::iterator itr = lower_bound( _testStats.begin(), _testStats.end(), _dataLLR );
  int nOver = distance( itr, _testStats.end() );

  return double( nOver ) / double( _testStats.size() );

}

void PValueTest::finalize() {

  double histMax = _testStats[_testStats.size() / 2] > _dataLLR ?
      2 * _testStats[_testStats.size() / 2] :
      2 * _dataLLR;
  if ( histMax < 1. ) histMax = 1.;
  int nBins = _testStats.size()/100;
  double off = 1.5 * histMax / double(nBins);

  string hName = str( format( "Likelihood_FCN-scale%2.2e" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(), "", nBins, -off, histMax-off );
  string xTitle = str( format( "-2ln#lambda( #Lambda = %2.2f TeV )" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );

  foreach( const double& x, _testStats )
    _minus2LnLikelihoodDistribution->Fill( x );
  TLine* dataLine = new TLine( _dataLLR, _minus2LnLikelihoodDistribution->GetMinimum(), _dataLLR,
                               _minus2LnLikelihoodDistribution->GetMaximum() );

  dataLine->SetLineColor( kRed );
  dataLine->SetLineWidth( 2 );

  TCanvas* pvc = new TCanvas( ( hName + "Canvas" ).c_str(), "", 500, 500 );
  pvc->cd();
  _minus2LnLikelihoodDistribution->Draw();
  dataLine->Draw();
  string cName = _minus2LnLikelihoodDistribution->GetName();
  pvc->Print( ( "figures/Likelihood/" + cName + ".pdf" ).c_str() );
}

double PValueTest::alpha() const {
  return _alpha;
}
void PValueTest::alpha( const double& alpha ) {
  _alpha = alpha;
}
