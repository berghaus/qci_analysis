#include "Likelihood.hpp"

#include <cfloat>
#include <stdexcept>
#include <TH1.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/CombinedMinimizer.h>
#include <Minuit2/MnStrategy.h>

#include "PDF.hpp"

using namespace std;
using namespace ROOT::Minuit2;

Likelihood_FCN::Likelihood_FCN() :
    _data( 0 ) {

  _pars.Add( "alpha", 0., 2.0e-2, 0., 4. );
}

Likelihood_FCN::Likelihood_FCN( const TH1* data, const PDF* pdf, const double alpha ) :
    _data( data ),
    _pdf( pdf ) {

  _pars.Add( "alpha", alpha, 2.0e-02, 0., 4. );

}

Likelihood_FCN::~Likelihood_FCN() {
}

double Likelihood_FCN::operator()() const {

  return ( *this )( _pars.Params() );

}

double Likelihood_FCN::operator()( const std::vector< double >& par ) const {
  const PDF& pdf = *_pdf;

  double result = 0.;
  for( int bin = 1; bin <= _data->GetNbinsX(); ++bin ) {
    double x = _data->GetBinCenter( bin );
    int n = _data->GetBinContent( bin );
    double logProb = pdf( x, n, par );
    result += -2 * logProb;
  }

  return result;
}

double Likelihood_FCN::Up() const {
  return 1.;
}

void Likelihood_FCN::data( const TH1* data ) {
  _data = data;
  _isMinimized = false;
}

void Likelihood_FCN::pdf( const PDF* pdf ) {
  _pdf = pdf;
  _isMinimized = false;
}

const TH1* Likelihood_FCN::data() const {
  return _data;
}
const PDF* Likelihood_FCN::pdf() const {
  return _pdf;
}
bool Likelihood_FCN::isMinimized() const {
  return _isMinimized;
}

double Likelihood_FCN::Minimize() {
  return Minimize( _pars );
}

double Likelihood_FCN::Minimize( MnUserParameters& pars ) {

  CombinedMinimizer combMin;
  MnStrategy strat( 2 );
  FunctionMinimum minResults( combMin.Minimize( *this, pars, strat, 10000, 0.1 ) );

  _pars = pars = minResults.UserParameters();
  _isMinimized = minResults.IsValid();
  return pars.Params().at( 0 );

}

vector< double > Likelihood_FCN::pars() {
  return _pars.Params();
}

void Likelihood_FCN::accept( TestStatMonitor& mon ) {
  mon.monitor( *this );
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

LikelihoodRatio_FCN::LikelihoodRatio_FCN() :
    _data( 0 ) {
}

LikelihoodRatio_FCN::LikelihoodRatio_FCN( const TH1* data, const PDF* pdf ) :
    _data( (TH1*) data->Clone( "data" ) ),
    _pdf( pdf ),
    _numerator( data, pdf ),
    _denominator( data, pdf ) {

  init();

}

LikelihoodRatio_FCN::~LikelihoodRatio_FCN() {
  _data->Delete();
}

void LikelihoodRatio_FCN::init() {
  _denominator.Minimize();
}

double LikelihoodRatio_FCN::operator()( const std::vector< double >& par ) {

  // would nomlize denominator over nuisance parameters .. but none for now
  if ( !_denominator.isMinimized() ) _denominator.Minimize();

  if ( _numerator( par ) < _denominator() ) {
    cout << "-2lnL(" << par.at( 0 ) << ") = " << _numerator( par ) << '\n' << "-2lnL(" << _denominator.pars().at( 0 )
         << ") = " << _denominator() << '\n';
    //throw( logic_error("minimized likelihood not at minimum") );

  }
  return _numerator( par ) - _denominator();

}

double LikelihoodRatio_FCN::Up() const {
  return 1.;
}

void LikelihoodRatio_FCN::data( const TH1* data ) {
  _data = (TH1*) data->Clone( "data" );
  _numerator.data( data );
  _denominator.data( data );
  init();
}

void LikelihoodRatio_FCN::pdf( const PDF* pdf ) {
  _pdf = pdf;
  _numerator.pdf( pdf );
  _denominator.pdf( pdf );
  init();
}

TH1* LikelihoodRatio_FCN::data() const {
  return _data;
}
const PDF* LikelihoodRatio_FCN::pdf() const {
  return _pdf;
}

void LikelihoodRatio_FCN::accept( TestStatMonitor& mon ) {
  mon.monitor( *this );
}
