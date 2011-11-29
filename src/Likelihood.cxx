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

#include "PseudoExperiment.hpp"
#include "PDF.hpp"

using namespace std;
using namespace ROOT::Minuit2;

Likelihood_FCN::Likelihood_FCN() :
    _data( 0 ) {

  _pars.Add( "alpha", 0., 1.e-3, 0., 16. );
}

Likelihood_FCN::Likelihood_FCN( const Experiment* data, const PDF* pdf, const double alpha ) :
    _data( data ),
    _pdf( pdf ) {

  _pars.Add( "alpha", alpha, 1.e-3, 0., 16. );

}

Likelihood_FCN::~Likelihood_FCN() {
}

double Likelihood_FCN::operator()() const {

  return ( *this )( _pars.Params() );

}

double Likelihood_FCN::operator()( const std::vector< double >& par ) const {

  const PDF& pdf = *_pdf;

  double result = 0.;
  for( int bin = 0; bin < _data->x().size(); ++bin ) {
    double x = _data->x( bin );
    int n = _data->y( bin );
    double logProb = pdf( x, n, par );
    result += -2 * logProb;
  }

  return result;
}

double Likelihood_FCN::Up() const {
  return 1.;
}

void Likelihood_FCN::data( const Experiment* data ) {
  _data = data;
  _isMinimized = false;
}

void Likelihood_FCN::pdf( const PDF* pdf ) {
  _pdf = pdf;
  _isMinimized = false;
}

const Experiment* Likelihood_FCN::data() const {
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
  if ( !_isMinimized ) cout << "invalid fit!\n";
  return 1. / pow( pars.Params().at( 0 ), 0.25 );

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

LikelihoodRatio::LikelihoodRatio() :
    _data( 0 ) {
}

LikelihoodRatio::LikelihoodRatio( const Experiment* data, const PDF* pdf, const double& scale ) :
  _data( data ),
    _pdf( pdf ),
    _numerator( data, pdf, scale ),
    _denominator( data, pdf, scale ) {

  init();

}

LikelihoodRatio::~LikelihoodRatio() {
  cout << "deleting a likelihoodratio\n";
  delete _pdf;
}

void LikelihoodRatio::init() {
  _denominator.Minimize();
}

double LikelihoodRatio::operator()( const std::vector< double >& par ) {

  // would nomlize denominator over nuisance parameters .. but none for now
  if ( !_denominator.isMinimized() ) _denominator.Minimize();
  return _numerator( par ) - _denominator();

}

double LikelihoodRatio::Up() const {
  return 1.;
}

void LikelihoodRatio::data( const Experiment* data ) {
  _data = data;
  _numerator.data( data );
  _denominator.data( data );
  init();
}

void LikelihoodRatio::pdf( const PDF* pdf ) {
  _pdf = pdf;
  _numerator.pdf( pdf );
  _denominator.pdf( pdf );
  init();
}

const Experiment* LikelihoodRatio::data() const {
  return _data;
}

const PDF* LikelihoodRatio::pdf() const {
  return _pdf;
}

void LikelihoodRatio::accept( TestStatMonitor& mon ) {
  mon.monitor( *this );
}

Likelihood_FCN LikelihoodRatio::numerator() const {
  return _numerator;
}

Likelihood_FCN LikelihoodRatio::denominator() const {
  return _denominator;
}
