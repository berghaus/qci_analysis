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
    _randomCompScale( 0.5, 20. ),
    _randomAlpha( 1., exp( pow( 1. / 3., 4 ) ) ) {
  init();
}

TestStatMonitor::TestStatMonitor( const string& folder, const string& ext ) :
    _folder( folder ),
    _ext( ext ),
    _randomCompScale( 0.5, 20. ),
    _randomAlpha( 1., exp( pow( 1. / 3., 4 ) ) ) {
  init();
}

void TestStatMonitor::init() {

  _minScale = 4.;
  _maxScale = 8.;
  _nBinsScale = 100;
  _randomCompScale( _minScale, _maxScale );

  double min = 0.1;
  double max = pow( 1. / 3., 4 ) + 0.1;
  double nBins = 1000;
  vector< double > alphaBins;
  alphaBins.reserve( nBins );
  for( int iBin = 0; iBin <= nBins; ++iBin ) {
    double binEdge = pow( max, double( iBin ) / double( nBins ) )
                     * pow( min, double( nBins - iBin ) / double( nBins ) );
    alphaBins.push_back( binEdge );
  }

  _likelihoodVsScale = new TH2D( "Likelihood", "likelihood", _nBinsScale, _minScale, _maxScale, 1000, -5., 200. );
  _likelihoodVsScale->SetXTitle( "#Lambda [TeV]" );
  _likelihoodVsScale->SetYTitle( "-2*ln( L(data|#Lambda) )" );

  _likelihoodRatioVsScale = new TH2D( "LikelihoodRatio", "likelihoodRatio", _nBinsScale, _minScale, _maxScale, 1000,
                                      -5., 200. );
  _likelihoodRatioVsScale->SetXTitle( "#Lambda [TeV]" );
  _likelihoodRatioVsScale->SetYTitle( "-2*ln( #lambda(#Lambda) )" );

  _likelihoodVsAlpha = new TH2D( "LikelihoodVsAlpha", "likelihood", alphaBins.size() - 1, &alphaBins[0], 1000, -5.,
                                 200. );
  _likelihoodVsAlpha->SetXTitle( "#alpha + 0.1 = #Lambda^{-4} + 0.1 [TeV^{-4}]" );
  _likelihoodVsAlpha->SetYTitle( "-2*ln( L(data|#alpha) )" );

  _likelihoodRatioVsAlpha = new TH2D( "LikelihoodRatioVsAlpha", "likelihood", alphaBins.size() - 1, &alphaBins[0], 1000,
                                      -5., 200. );
  _likelihoodRatioVsAlpha->SetXTitle( "#alpha + 0.1 = #Lambda^{-4} + 0.1 [TeV^{-4}]" );
  _likelihoodRatioVsAlpha->SetYTitle( "-2*ln( #lambda(#alpha) )" );

  _minimizedAlpha = new TH1D( "MinimizedAlpha", "minimizedAlpha", 1700, -0.5, 16.5 );
  _minimizedAlpha->SetXTitle( "#alpha = 1/#Lambda^{4} [TeV^{-4}]" );
  _minimizedAlpha->SetYTitle( "Number of PEs" );

  _minimizedLaunda = new TH1D( "MinimizedLaunda", "minimizedLaunda", _nBinsScale, _minScale, _maxScale );
  _minimizedLaunda->SetXTitle( "#Lambda [TeV]" );
  _minimizedLaunda->SetYTitle( "Number of PEs" );

}

void TestStatMonitor::finalize() {

  TCanvas *lc = new TCanvas( "LikelihoodCanvas", "", 500, 500 );
  lc->cd();
  //lc->SetLogy();
  _likelihoodVsScale->Draw( "COLZ" );
  lc->Print( ( _folder + string( _likelihoodVsScale->GetName() ) + _ext ).c_str() );

  TCanvas *lrc = new TCanvas( "LikelihoodRatioCanvas", "", 500, 500 );
  lrc->cd();
  //lrc->SetLogy();
  _likelihoodRatioVsScale->Draw( "COLZ" );
  lrc->Print( ( _folder + string( _likelihoodRatioVsScale->GetName() ) + _ext ).c_str() );

  TCanvas *likelihoodVsAlphaCanvas = new TCanvas( "LikelihoodVsAlphaCanvas", "", 500, 500 );
  likelihoodVsAlphaCanvas->cd();
  likelihoodVsAlphaCanvas->SetLogx();
  _likelihoodVsAlpha->Draw( "COLZ" );
  likelihoodVsAlphaCanvas->Print( ( _folder + string( _likelihoodVsAlpha->GetName() ) + _ext ).c_str() );

  TCanvas *likelihoodRatioVsAlphaCanvas = new TCanvas( "LikelihoodRatioVsAlphaCanvas", "", 500, 500 );
  likelihoodRatioVsAlphaCanvas->cd();
  likelihoodRatioVsAlphaCanvas->SetLogx();
  _likelihoodRatioVsAlpha->Draw( "COLZ" );
  likelihoodRatioVsAlphaCanvas->Print( ( _folder + string( _likelihoodRatioVsAlpha->GetName() ) + _ext ).c_str() );

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

  _likelihoodVsScale->Delete();
  _likelihoodRatioVsScale->Delete();
  _minimizedAlpha->Delete();
  _minimizedLaunda->Delete();
}

void TestStatMonitor::monitor( Likelihood_FCN& l ) {

  for( int i = 0; i < 100; ++i ) {
    double scale = _randomCompScale();
    double alpha = log( _randomAlpha() );
    _likelihoodVsScale->Fill( scale, l( vector< double >( 1, 1. / pow( scale, 4 ) ) ) );
    _likelihoodVsAlpha->Fill( alpha + 0.1, l( vector< double >( 1, alpha ) ) );
  }

  if ( l.isMinimized() ) {
    _minimizedAlpha->Fill( l.pars().at( 0 ) );
    _minimizedLaunda->Fill( pow( 1. / l.pars().at( 0 ), 0.25 ) );
  }

}

void TestStatMonitor::monitor( LikelihoodRatio& launda ) {

  for( int i = 0; i < 100; ++i ) {
    double scale = _randomCompScale();
    double alpha = log( _randomAlpha() );
    _likelihoodRatioVsScale->Fill( scale, launda( vector< double >( 1, 1. / pow( scale, 4 ) ) ) );
    _likelihoodRatioVsAlpha->Fill( alpha + 0.1, launda( vector< double >( 1, alpha ) ) );
  }

}

