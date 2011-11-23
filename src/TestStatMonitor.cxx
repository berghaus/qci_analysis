#include "TestStatMonitor.hpp"
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

#include "Likelihood.hpp"

using namespace std;
using namespace boost::assign;

TestStatMonitor::TestStatMonitor() :
    _folder( "figures/" ),
    _ext( ".png" ),
    _randomCompScale( 0.5, 20. ) {
  init();
}

TestStatMonitor::TestStatMonitor( const string& folder, const string& ext ) :
    _folder( folder ),
    _ext( ext ),
    _randomCompScale( 0.5, 20. ) {
  init();
}

void TestStatMonitor::init() {

  _minScale = 4.;
  _maxScale = 8.;
  _nBinsScale = 100;
  _randomCompScale( _minScale, _maxScale );

  double min = 40;
  double max = 1.e3;
  double nBins = 1000;
  vector< double > logLBins;
  logLBins.reserve( nBins );
  for( int iBin = 0; iBin <= nBins; ++iBin ) {
    double binEdge = pow( max, double( iBin ) / double( nBins ) )
                     * pow( min, double( nBins - iBin ) / double( nBins ) );
    logLBins.push_back( binEdge );
  }

  min = 0.01;
  vector< double > logLambdaBins;
  for( int iBin = 0; iBin <= nBins; ++iBin ) {
    double binEdge = pow( max, double( iBin ) / double( nBins ) )
                     * pow( min, double( nBins - iBin ) / double( nBins ) );
    logLambdaBins.push_back( binEdge );
  }

  _likelihood = new TH2D( "Likelihood", "likelihood", _nBinsScale, _minScale, _maxScale, logLBins.size() - 1, &logLBins[0] );
  _likelihood->SetXTitle( "#Lambda [TeV]" );
  _likelihood->SetYTitle( "-2*ln( L(data|#Lambda) )" );

  _likelihoodRatio = new TH2D( "LikelihoodRatio", "likelihoodRatio", _nBinsScale, _minScale, _maxScale, logLambdaBins.size() - 1,
                               &logLambdaBins[0] );
  _likelihoodRatio->SetXTitle( "#Lambda [TeV]" );
  _likelihoodRatio->SetYTitle( "-2*ln( #lambda(#Lambda) )" );

  _minimizedAlpha = new TH1D( "MinimizedAlpha", "minimizedAlpha", 200, 0., -1. );
  _minimizedAlpha->SetXTitle( "#alpha = 1/#Lambda^{4} [TeV^{-4}]" );
  _minimizedAlpha->SetYTitle( "Number of PEs" );

  _minimizedLaunda = new TH1D( "MinimizedLaunda", "minimizedLaunda", _nBinsScale, _minScale, _maxScale );
  _minimizedLaunda->SetXTitle( "#Lambda [TeV]" );
  _minimizedLaunda->SetYTitle( "Number of PEs" );

}

void TestStatMonitor::finalize() {

  TCanvas *lc = new TCanvas( "LikelihoodCanvas", "", 500, 500 );
  lc->cd();
  lc->SetLogy();
  _likelihood->Draw( "COLZ" );
  lc->Print( ( _folder + string( _likelihood->GetName() ) + _ext ).c_str() );

  TCanvas *lrc = new TCanvas( "LikelihoodRatioCanvas", "", 500, 500 );
  lrc->cd();
  lrc->SetLogy();
  _likelihoodRatio->Draw( "COLZ" );
  lrc->Print( ( _folder + string( _likelihoodRatio->GetName() ) + _ext ).c_str() );

  TCanvas *ac = new TCanvas( "AlphaCanvas", "", 500, 500 );
  ac->cd();
  _minimizedAlpha->Draw();
  ac->Print( ( _folder + string( _minimizedAlpha->GetName() ) + _ext ).c_str() );

  TCanvas *sc = new TCanvas( "LaundaCanvas", "", 500, 500 );
  sc->cd();
  _minimizedLaunda->Draw();
  sc->Print( ( _folder + string( _minimizedLaunda->GetName() ) + _ext ).c_str() );

}

TestStatMonitor::~TestStatMonitor() {

  _likelihood->Delete();
  _likelihoodRatio->Delete();
  _minimizedAlpha->Delete();
  _minimizedLaunda->Delete();
}

void TestStatMonitor::monitor( Likelihood_FCN& l ) {

  for( int i = 0; i < 1000; ++i ) {
    double scale = _randomCompScale();
    double alpha = 1. / pow( scale, 4 );
    _likelihood->Fill( scale, l( vector< double >( 1, alpha ) ) );
  }

  if ( l.isMinimized() ) {
    _minimizedAlpha->Fill( l.pars().at( 0 ) );
    _minimizedLaunda->Fill( 1. / pow( l.pars().at( 0 ), 0.25 ) );
  }

}

void TestStatMonitor::monitor( LikelihoodRatio& launda ) {

  for( int i = 0; i < 100; ++i ) {
    double scale = _randomCompScale();
    double alpha = 1. / pow( scale, 4 );
    _likelihoodRatio->Fill( scale, launda( vector< double >( 1, alpha ) ) );
  }

}

