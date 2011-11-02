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

TestStatMonitor::TestStatMonitor()
  : _folder("figures/")
  , _ext   (".png") {
  init();
}


TestStatMonitor::TestStatMonitor( const string& folder, const string& ext )
  : _folder( folder )
  , _ext   ( ext ) {
  init();
}


void
TestStatMonitor::init() {

  _likelihood      = new TH2D("Likelihood"     ,"likelihood",     500, 0.,4.e-6,500,0,-1);
  _likelihoodRatio = new TH2D("LikelihoodRatio","likelihoodRatio",500, 0.,4.e-6,500,0,-1);
  _minimizedAlpha  = new TH1D("MinimizedAlpha" ,"minimizedAlpha", 400, 0.,4.e-6 );

}


TestStatMonitor::~TestStatMonitor() {

  TCanvas c("TestStatCanvas","",800,600); c.cd();

  _likelihood     ->Draw("COLZ");
  c.Print( (_folder+string( _likelihood->GetName() )+_ext).c_str() );

  _likelihoodRatio->Draw("COLZ");
  c.Print( (_folder+string( _likelihoodRatio->GetName() )+_ext).c_str() );

  _minimizedAlpha ->Draw();
  c.Print( (_folder+string( _minimizedAlpha->GetName() )+_ext).c_str() );
  
  _likelihood     ->Delete();
  _likelihoodRatio->Delete();
  _minimizedAlpha ->Delete();
}


void
TestStatMonitor::monitor( Likelihood_FCN& l ) {

  for( int i = 0; i < 100; ++i ) {
    double randAlpha = 4.e-6 * double( rand() % 1000 ) / 1000.;
    _likelihood->Fill( randAlpha, l( vector<double>(1, randAlpha ) ) );
  }
  
  if ( l.isMinimized() )
    _minimizedAlpha->Fill( l.pars().at(0) );

}


void
TestStatMonitor::monitor( LikelihoodRatio_FCN& launda ) {
  for( int i = 0; i < 100; ++i ) {
    double randAlpha = 4.e-6 * double( rand() % 1000 ) / 1000.;
    _likelihoodRatio->Fill( randAlpha, launda( vector<double>(1, randAlpha ) ) );
  }

}

