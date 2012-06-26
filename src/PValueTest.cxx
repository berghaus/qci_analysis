#include "PValueTest.hpp"
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/checked_delete.hpp>

#include <TCanvas.h>
#include <TFitterMinuit.h>
#include <TH1.h>
#include <TLine.h>

#include "Prediction.hpp"
#include "Likelihood.hpp"

#define foreach BOOST_FOREACH
using namespace std;
using boost::format;

PValueTest::PValueTest( const double alpha, const vector< Neg2LogMaximumLikelihoodRatio* >& lambdas ) :
    _alpha( alpha ),
    _lambdas( lambdas ),
    _dataLLR( 0 ) {
  init();
}

PValueTest::PValueTest( const double& alpha, const int& nPE, PseudoExperimentFactory& peFactory ) :
    _alpha( alpha ),
    _dataLLR( 0 ) {
  init( nPE, peFactory );
}

void PValueTest::init( const int& nPE, PseudoExperimentFactory& peFactory ) {

  vector< double > par( 1, _alpha );
  for( int i = 0; i < nPE; ++i ) {
    PseudoExperiment pe = peFactory.build( _alpha );
    Prediction * pePDF = new Prediction( *peFactory.pdf() );
    pePDF->nData( pe );
    Neg2LogMaximumLikelihoodRatio n2llr( &pe, pePDF, _alpha );
    for( double scale = 2.; scale < 8.; scale += 0.1 )
      n2llr( vector< double >( 1, scale ) );

    _testStats.push_back( n2llr( par ) );
  }

}

void PValueTest::init() {

  vector< double > par( 1, _alpha );
  _testStats.reserve( _lambdas.size() );
  foreach( Neg2LogMaximumLikelihoodRatio* lambda, _lambdas )
  {
    Neg2LogMaximumLikelihoodRatio& l = *lambda;
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

double PValueTest::operator()( Neg2LogMaximumLikelihoodRatio& lambda ) {

  vector< double > par( 1, _alpha );

  return ( *this )( lambda( par ) );

}

double PValueTest::operator()( const double& llr ) {

  _dataLLR = llr;
  vector< double >::iterator itr = lower_bound( _testStats.begin(), _testStats.end(), _dataLLR );
  int nOver = distance( itr, _testStats.end() );

  return double( nOver ) / double( _testStats.size() );

}

void PValueTest::finalize( const std::string& dir ) {

  double histMax = _testStats[_testStats.size() / 2] > _dataLLR ?
      2 * _testStats[_testStats.size() / 2] :
      2 * _dataLLR;
  if ( histMax < 1. ) histMax = 1.;
  int nBins = _testStats.size() / 100;
  double off = 1.5 * histMax / double( nBins );

  string hName = str( format( "Likelihood_FCN-scale%2.2e" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(), "", nBins, -off, histMax - off );
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
  pvc->Print( ( dir + cName + ".eps" ).c_str() );
}

double PValueTest::alpha() const {
  return _alpha;
}
void PValueTest::alpha( const double& alpha ) {
  _alpha = alpha;
}

void PValueTest::clear() {

  _alpha = -1.;
  _lambdas.clear();
  _dataLLR = 0.;

  //delete _minus2LnLikelihoodDistribution;
  _testStats.clear();

}

ostream& operator<<( ostream& out, const PValueTest& test ) {

  out << test.alpha() << " ";
  vector< double >::const_iterator itr = test._testStats.begin();
  vector< double >::const_iterator end = test._testStats.end();
  for( ; itr != end; ++itr ) {
    out << *itr << " ";
  }

  return out;
}

istream& operator>>( istream& in, PValueTest& test ) {

  test.clear();

  string buffer;
  double x;

  getline( in, buffer );
  if ( !in ) return in;
  stringstream bufferStream( buffer );
  bufferStream >> x;
  test.alpha( x );
  while ( bufferStream >> x ) {
    test._testStats.push_back( x );
  }

  return in;
}
