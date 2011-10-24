#include "PValueTest.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>

#include <boost/format.hpp>
#include <boost/checked_delete.hpp>

#include <TFitterMinuit.h>
#include <TH1.h>

#include "PDF.hpp"
#include "Likelihood.hpp"

using namespace std;
using boost::format;


PValueTest::  PValueTest( const double alpha, const LikelihoodRatio& testStat, PseudoExperimentFactory* factory )
  : _alpha   ( alpha )
  , _testStat( testStat )
  , _factory ( factory ) {
}


void
PValueTest::init( TFitterMinuit* fitter, const PDF* pdf, unsigned int nPE ) {

  _nPE = nPE;

  string hName = str( format("Likelihood-alpha%2.1e") % _alpha  );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(),"",500,0.,-1.);
  string xTitle = str( format("-2ln#Lampda(%2.1e)") % _alpha  );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );

  PseudoExperiment* pe = 0;
  for( int i = 0; i < _nPE; ++i ) {
    pe = _factory->build( _alpha );
    LikelihoodRatio * launda = new LikelihoodRatio( fitter, pe, pdf );
    double min2LogLaunda = _testStat( pe );
    _minus2LnLikelihoodDistribution->Fill( min2LogLaunda );
    _pes.push_back( launda );
  }
}


PValueTest::PValueTest( )
  : _alpha( 0 )
  , _factory ( 0 ) {
}


PValueTest::~PValueTest() {
  for_each( _pes.begin(), _pes.end(), boost::checked_deleter<LikelihoodRatio>() );
}

double
PValueTest::operator() ( const TH1* data ) {

  if ( !_factory ) throw( runtime_error("must provide PseudoExperimentFactory") );

  

  return 0.;

}


