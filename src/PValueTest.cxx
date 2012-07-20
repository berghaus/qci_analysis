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


template< typename Likelihood >
PValueTest<Likelihood>::PValueTest( const double alpha, const vector< Likelihood* >& lambdas ) :
    _alpha( alpha ),
    _lambdas( lambdas ),
    _dataLLR( 0 ) {
  init();
}


template< typename Likelihood >
void PValueTest<Likelihood>::init() {

  vector< double > par( 1, _alpha );
  _testStats.reserve( _lambdas.size() );
  foreach( Likelihood* lambda, _lambdas )
  {
    Likelihood& l = *lambda;
    _testStats.push_back( l( par ) );
  }
  sort( _testStats.begin(), _testStats.end() );

}

template< typename Likelihood >
PValueTest<Likelihood>::PValueTest() :
    _alpha( 0 ),
    _dataLLR( 0 ) {
}

template< typename Likelihood >
PValueTest<Likelihood>::~PValueTest() {
}

template< typename Likelihood >
double PValueTest<Likelihood>::operator()( Likelihood& lambda ) {

  vector< double > par( 1, _alpha );

  return ( *this )( lambda( par ) );

}

template< typename Likelihood >
double PValueTest<Likelihood>::operator()( const double& llr ) {

  _dataLLR = llr;
  vector< double >::iterator itr = lower_bound( _testStats.begin(), _testStats.end(), _dataLLR );
  int nOver = distance( itr, _testStats.end() );

  return double( nOver ) / double( _testStats.size() );

}

template< typename Likelihood >
void PValueTest<Likelihood>::finalize( const std::string& dir ) {

  double histMax = _testStats[_testStats.size() / 2] > _dataLLR ?
      2 * _testStats[_testStats.size() / 2] :
      2 * _dataLLR;
  if ( histMax < 1. ) histMax = 1.;
  int nBins = _testStats.size() / 100;
  double off = 1.5 * histMax / double( nBins );

  string hName = str( format( "Likelihood_FCN-scale%2.2e" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(), "", nBins, -off, histMax - off );
  string xTitle = _alpha == 0 ? "-2ln#lambda( #Lambda = #infty TeV )"
                              : str( format( "-2ln#lambda( #Lambda = %2.2f TeV )" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );

  foreach( const double& x, _testStats )
    _minus2LnLikelihoodDistribution->Fill( x );
  TLine* dataLine = new TLine( _dataLLR, _minus2LnLikelihoodDistribution->GetMinimum(), _dataLLR,
                               _minus2LnLikelihoodDistribution->GetMaximum() * 2 );

  dataLine->SetLineColor( kRed );
  dataLine->SetLineWidth( 2 );

  TCanvas* pvc = new TCanvas( ( hName + "Canvas" ).c_str(), "", 500, 500 );
  pvc->SetLogy();
  pvc->cd();
  _minus2LnLikelihoodDistribution->Draw();
  dataLine->Draw();
  string cName = _minus2LnLikelihoodDistribution->GetName();
  pvc->Print( ( dir + cName + ".pdf" ).c_str() );
}

template< typename Likelihood >
double PValueTest<Likelihood>::alpha() const {
  return _alpha;
}

template< typename Likelihood >
void PValueTest<Likelihood>::alpha( const double& alpha ) {
  _alpha = alpha;
}

template< typename Likelihood >
void PValueTest<Likelihood>::clear() {

  _alpha = -1.;
  _lambdas.clear();
  _dataLLR = 0.;

  //delete _minus2LnLikelihoodDistribution;
  _testStats.clear();

}

template< typename Likelihood >
ostream& operator<<( ostream& out, const PValueTest<Likelihood>& test ) {

  out << test.alpha() << " ";
  vector< double >::const_iterator itr = test._testStats.begin();
  vector< double >::const_iterator end = test._testStats.end();
  for( ; itr != end; ++itr ) {
    out << *itr << " ";
  }

  return out;
}


template< typename Likelihood >
istream& operator>>( istream& in, PValueTest<Likelihood>& test ) {

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


template class PValueTest<Neg2LogLikelihood_FCN>;
template typename std::ostream& operator<< <Neg2LogLikelihood_FCN>( std::ostream&, const PValueTest<Neg2LogLikelihood_FCN>& );
template typename std::istream& operator>> <Neg2LogLikelihood_FCN>( std::istream&, PValueTest<Neg2LogLikelihood_FCN>& );

template class PValueTest<Neg2LogMaximumLikelihoodRatio>;
template typename std::ostream& operator<< <Neg2LogMaximumLikelihoodRatio>( std::ostream&, const PValueTest<Neg2LogMaximumLikelihoodRatio>& );
template typename std::istream& operator>> <Neg2LogMaximumLikelihoodRatio>( std::istream&, PValueTest<Neg2LogMaximumLikelihoodRatio>& );

