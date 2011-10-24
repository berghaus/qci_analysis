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


PValueTest::PValueTest( const double alpha, const LikelihoodRatio& testStat, PseudoExperimentFactory* factory )
  : _alpha   ( alpha )
  , _testStat( testStat )
  , _factory ( factory ) {
  init();
}


void
PValueTest::init( unsigned int nPE ) {

  _nPE = nPE;

  string hName = str( format("Likelihood-alpha%2.1e") % _alpha  );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(),"",100,0.,-1.);
  string xTitle = str( format("-2ln#lambda(%2.1e GeV^{-2})") % _alpha  );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );

  vector<double> par(1,_alpha);
  PseudoExperiment* pe = 0;
  for( int i = 0; i < _nPE; ++i ) {

    pe = _factory->build( _alpha );
    _testStat.data( pe );
    _minus2LnLikelihoodDistribution->Fill( _testStat( par ) );
  }
}


PValueTest::PValueTest( )
  : _alpha( 0 )
  , _factory ( 0 ) {
}


PValueTest::~PValueTest() {
  //for_each( _pes.begin(), _pes.end(), boost::checked_deleter<LikelihoodRatio>() );
  TCanvas c("c","",800,600);
  _minus2LnLikelihoodDistribution->Draw();
  string cName = _minus2LnLikelihoodDistribution->GetName();
  c.Print( (cName+".png").c_str() );
}

double
PValueTest::operator() ( const TH1* data ) {

  if ( !_factory ) throw( runtime_error("must provide PseudoExperimentFactory") );

  vector<double> par( 1, _alpha );
  _testStat.data( data );
  int dataOutcomeBin = _minus2LnLikelihoodDistribution->FindBin( _testStat( par ) );

  return _minus2LnLikelihoodDistribution->Integral(dataOutcomeBin,-1) / _minus2LnLikelihoodDistribution->Integral();

}


