#include "PValueTest.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>

#include <boost/format.hpp>
#include <boost/checked_delete.hpp>

#include <TCanvas.h>
#include <TFitterMinuit.h>
#include <TH1.h>

#include "PDF.hpp"
#include "Likelihood.hpp"

using namespace std;
using boost::format;

PValueTest::PValueTest( const double alpha, const LikelihoodRatio_FCN& testStat, vector<PseudoExperiment*> pes)
  : _alpha   ( alpha )
  , _testStat( testStat )
  , _pes ( pes )
  , _testStats( pes.size() ) {
  init();
}


void
PValueTest::init() {


  string hName = str( format("Likelihood_FCN-alpha%2.1e") % _alpha  );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(),"",1000,0.,-1.);
  string xTitle = str( format("-2ln#lambda(%2.1e GeV^{-2})") % _alpha  );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );

  vector<double> par(1,_alpha);
  vector<PseudoExperiment*>::iterator itr = _pes.begin();
  vector<PseudoExperiment*>::iterator end = _pes.end();
  for( ; itr != end; ++itr ) {
    PseudoExperiment* pe = *itr;
    _testStat.data( pe );
    _testStats.push_back( _testStat( par ) );
    _minus2LnLikelihoodDistribution->Fill( _testStat( par ) );
  }
  sort( _testStats.begin(), _testStats.end() );
}


PValueTest::PValueTest( )
  : _alpha( 0 ) {
}


PValueTest::~PValueTest() {
  //for_each( _pes.begin(), _pes.end(), boost::checked_deleter<LikelihoodRatio_FCN>() );
  TCanvas* c = new TCanvas("PcalueCanvas","",500,500); c->cd();
  _minus2LnLikelihoodDistribution->Draw();
  string cName = _minus2LnLikelihoodDistribution->GetName();
  c->Print( (cName+".png").c_str() );
}

double
PValueTest::operator() ( const TH1* data ) {

  vector<double> par( 1, _alpha );
  _testStat.data( data );
  double dataLL = _testStat( par );
  vector<double>::iterator itr = lower_bound( _testStats.begin(), _testStats.end(), dataLL );
  int nOver = distance( itr, _testStats.end() );

  return  double(nOver) / double(_testStats.size());

}


