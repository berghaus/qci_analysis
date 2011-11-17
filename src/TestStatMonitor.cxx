#include "TestStatMonitor.hpp"
#include <cmath>
#include <cstdlib>
#include <string>

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
    _ext( ".png" ) {
  init();
}

TestStatMonitor::TestStatMonitor( const string& folder, const string& ext ) :
    _folder( folder ),
    _ext( ext ) {
  init();
}

void TestStatMonitor::finalize() {

  TCanvas *lc = new TCanvas( "LikelihoodCanvas", "", 500, 500 );
  lc->cd();
  _likelihood->Draw( "COLZ" );
  lc->Print( ( _folder + string( _likelihood->GetName() ) + _ext ).c_str() );

  TCanvas *lrc = new TCanvas( "LikelihoodRatioCanvas", "", 500, 500 );
  lrc->cd();
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

void TestStatMonitor::init() {

  _likelihood = new TH2D( "Likelihood", "likelihood", 2000, 0.4, 20., 5000, 0, -1 );
  _likelihood->SetXTitle( "#Lambda [TeV]" );
  _likelihood->SetYTitle( "L(data|#Lambda)" );

  _likelihoodRatio = new TH2D( "LikelihoodRatio", "likelihoodRatio", 2000, 0.4, 20., 5000, 0, -1 );
  _likelihoodRatio->SetXTitle( "#Lambda [TeV]" );
  _likelihoodRatio->SetYTitle( "#lambda(#Lambda)" );

  _minimizedAlpha = new TH1D( "MinimizedAlpha", "minimizedAlpha", 2000, 0., 16. );
  _minimizedAlpha->SetXTitle( "#alpha = 1/#Lambda^{4} [TeV^{-4}]" );
  _minimizedAlpha->SetYTitle( "Number of PEs" );

  _minimizedLaunda = new TH1D( "MinimizedLaunda", "minimizedLaunda", 2000, 0., -1 );
  _minimizedLaunda->SetXTitle( "#Lambda [TeV]" );
  _minimizedLaunda->SetYTitle( "Number of PEs" );

}

TestStatMonitor::~TestStatMonitor() {

  _likelihood->Delete();
  _likelihoodRatio->Delete();
  _minimizedAlpha->Delete();
  _minimizedLaunda->Delete();
}

void TestStatMonitor::monitor( Likelihood_FCN& l ) {

  for( int i = 0; i < 1000; ++i ) {
    double randAlpha = 4. * double( rand() % 2000 ) / 2000.;
    _likelihood->Fill( randAlpha, l( vector< double >( 1, randAlpha ) ) );
  }

  if ( l.isMinimized() ) {
    _minimizedAlpha->Fill( l.pars().at( 0 ) );
    _minimizedLaunda->Fill( 1./pow( l.pars().at( 0 ), 0.25 ) );
  }

}

void TestStatMonitor::monitor( LikelihoodRatio_FCN& launda ) {
  for( int i = 0; i < 100; ++i ) {
    double randAlpha = 4. * double( rand() % 2000 ) / 2000.;
    _likelihoodRatio->Fill( randAlpha, launda( vector< double >( 1, randAlpha ) ) );
  }

}

