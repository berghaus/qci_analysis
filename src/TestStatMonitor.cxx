#include "TestStatMonitor.hpp"
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

#include "Likelihood.hpp"

using namespace std;
using boost::format;
using namespace boost::assign;

TestStatMonitor::TestStatMonitor() :
    _alpha( 0 ),
    _label( "" ),
    _folder( "figures/" ),
    _ext( ".pdf" ),
    _randomCompScale( 0.5, 20. ),
    _randomAlpha( log( 0.001 ), log( 16. ) ) {
  init();
}

TestStatMonitor::TestStatMonitor( const double& alpha, const string& folder, const string& ext ) :
    _alpha( alpha ),
    _folder( folder ),
    _ext( ext ),
    _randomCompScale( 5., 20. ),
    _randomAlpha( log( 0.001 ), log( 16. ) ) {
  _label = _alpha < 0. ?
      "data" :
      ( str( format( "scale%2.1f" ) % pow( alpha, -0.25 ) ) );
  init();
}

void TestStatMonitor::init() {

  if ( _alpha <= 0. ) {
    _minScale = 5.;
    _maxScale = 20.;
    _nBinsScale = 150;
  } else {
    _minScale = pow( _alpha, -0.25 ) / 2;
    _maxScale = pow( _alpha, -0.25 ) * 2;
    _nBinsScale = 100;
  }
  _randomCompScale( _minScale, _maxScale );

  double min = 0.1;
  double max = pow( 1. / .5, 4 ) + min;
  double nBins = 2000;
  vector< double > alphaBins;
  alphaBins.reserve( nBins );
  for( int iBin = 0; iBin <= nBins; ++iBin ) {
    double binEdge = pow( max, double( iBin ) / double( nBins ) )
                     * pow( min, double( nBins - iBin ) / double( nBins ) );
    alphaBins.push_back( binEdge - min );
  }

  _likelihoodVsScale = new TH2D( ( _label + "Likelihood" ).c_str(), "likelihood", _nBinsScale, _minScale, _maxScale,
                                 1000, -10., 200. );
  _likelihoodVsScale->SetXTitle( "#Lambda [TeV]" );
  _likelihoodVsScale->SetYTitle( ( "-2*ln( L(" + _label + "|#Lambda) )" ).c_str() );

  _likelihoodRatioVsScale = new TH2D( ( _label + "LikelihoodRatio" ).c_str(), "likelihoodRatio", _nBinsScale, _minScale,
                                      _maxScale, 1000, -10., 200. );
  _likelihoodRatioVsScale->SetXTitle( "#Lambda [TeV]" );
  _likelihoodRatioVsScale->SetYTitle(
      str( format( "-2*ln( #lambda(#Lambda = %2.1f ) )" ) % pow( _alpha, -0.25 ) ).c_str() );

  _likelihoodVsAlpha = new TH2D( ( _label + "LikelihoodVsAlpha" ).c_str(), "likelihood", alphaBins.size() - 1,
                                 &alphaBins[0], 1000, -10., 200. );
  _likelihoodVsAlpha->SetXTitle( "#alpha = #Lambda^{-4} [TeV^{-4}]" );
  _likelihoodVsAlpha->SetYTitle( ( "-2*ln( L(" + _label + "|#alpha) )" ).c_str() );

  _likelihoodRatioVsAlpha = new TH2D( ( _label + "LikelihoodRatioVsAlpha" ).c_str(), "likelihood", alphaBins.size() - 1,
                                      &alphaBins[0], 1000, -10., 200. );
  _likelihoodRatioVsAlpha->SetXTitle( "#alpha = #Lambda^{-4} [TeV^{-4}]" );
  _likelihoodRatioVsAlpha->SetYTitle(
      str( format( "-2*ln( #lambda(#Lambda = %2.1f ) )" ) % pow( _alpha, -0.25 ) ).c_str() );

  _minimizedAlpha = new TH1D( ( _label + "MinimizedAlpha" ).c_str(), "minimizedAlpha", 1000., 0., -1. );
  _minimizedAlpha->SetXTitle( "#alpha = #Lambda^{-4} [TeV^{-4}]" );
  _minimizedAlpha->SetYTitle( "Number of PEs" );

  _minimizedLaunda = new TH1D( ( _label + "MinimizedLaunda" ).c_str(), "minimizedLaunda", 1000., 0., -1. );
  _minimizedLaunda->SetXTitle( "#Lambda [TeV]" );
  _minimizedLaunda->SetYTitle( "Number of PEs" );

}

void TestStatMonitor::finalize() {

  TCanvas *lc = new TCanvas( ( _label + "LikelihoodCanvas" ).c_str(), "", 500, 500 );
  lc->cd();
  //lc->SetLogy();
  _likelihoodVsScale->Draw( "COLZ" );
  lc->Print( ( _folder + string( _likelihoodVsScale->GetName() ) + _ext ).c_str() );

  TCanvas *lrc = new TCanvas( ( _label + "LikelihoodRatioCanvas" ).c_str(), "", 500, 500 );
  lrc->cd();
  //lrc->SetLogy();
  _likelihoodRatioVsScale->Draw( "COLZ" );
  lrc->Print( ( _folder + string( _likelihoodRatioVsScale->GetName() ) + _ext ).c_str() );

  TCanvas *likelihoodVsAlphaCanvas = new TCanvas( ( _label + "LikelihoodVsAlphaCanvas" ).c_str(), "", 500, 500 );
  likelihoodVsAlphaCanvas->cd();
  likelihoodVsAlphaCanvas->SetLogx();
  _likelihoodVsAlpha->Draw( "COLZ" );
  likelihoodVsAlphaCanvas->Print( ( _folder + string( _likelihoodVsAlpha->GetName() ) + _ext ).c_str() );

  TCanvas *likelihoodRatioVsAlphaCanvas = new TCanvas( ( _label + "LikelihoodRatioVsAlphaCanvas" ).c_str(), "", 500,
                                                       500 );
  likelihoodRatioVsAlphaCanvas->cd();
  likelihoodRatioVsAlphaCanvas->SetLogx();
  _likelihoodRatioVsAlpha->Draw( "COLZ" );
  likelihoodRatioVsAlphaCanvas->Print( ( _folder + string( _likelihoodRatioVsAlpha->GetName() ) + _ext ).c_str() );

  TCanvas *ac = new TCanvas( ( _label + "AlphaCanvas" ).c_str(), "", 500, 500 );
  ac->cd();
  ac->SetLogx();
  _minimizedAlpha->Draw();
  ac->Print( ( _folder + string( _minimizedAlpha->GetName() ) + _ext ).c_str() );

  TCanvas *sc = new TCanvas( ( _label + "LaundaCanvas" ).c_str(), "", 500, 500 );
  sc->cd();
  _minimizedLaunda->Draw();
  sc->Print( ( _folder + string( _minimizedLaunda->GetName() ) + _ext ).c_str() );

}

TestStatMonitor::~TestStatMonitor() {

}

void TestStatMonitor::monitor( Neg2LogLikelihood_FCN& l ) {

  for( int i = 0; i < 1000; ++i ) {
    double scale = _randomCompScale();
    double alpha = exp( _randomAlpha() ) - 0.001;
    _likelihoodVsScale->Fill( scale, l( vector< double >( 1, 1. / pow( scale, 4 ) ) ) );
    _likelihoodVsAlpha->Fill( alpha, l( vector< double >( 1, alpha ) ) );
  }

  if ( l.isMinimized() ) {
    _minimizedAlpha->Fill( l.pars().at( 0 ) );
    _minimizedLaunda->Fill( pow( 1. / l.pars().at( 0 ), 0.25 ) );
  }

}

void TestStatMonitor::monitor( Neg2LogLikelihoodRatio& launda ) {

  for( int i = 0; i < 1000; ++i ) {
    double scale = _randomCompScale();
    double alpha = exp( _randomAlpha() ) - 0.001;
    _likelihoodRatioVsScale->Fill( scale, launda( vector< double >( 1, 1. / pow( scale, 4 ) ) ) );
    _likelihoodRatioVsAlpha->Fill( alpha, launda( vector< double >( 1, alpha ) ) );
  }

}

