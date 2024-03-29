#include "Likelihood.hpp"

#include <cfloat>
#include <stdexcept>
#include <boost/foreach.hpp>
#ifndef foreach
#define foreach BOOST_FOREACH
#endif
#include <boost/lexical_cast.hpp>
#include <TH1.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/CombinedMinimizer.h>
#include <Minuit2/MnStrategy.h>

#include "PseudoExperiment.hpp"
#include "Prediction.hpp"
#include "TestStatMonitor.hpp"

using namespace std;
using namespace ROOT::Minuit2;
using boost::lexical_cast;

Neg2LogLikelihood_FCN::Neg2LogLikelihood_FCN() :
  _data( 0 ) {

  _pars.Add( "alpha", 0., 2.e-4, 0., 16. );
}

Neg2LogLikelihood_FCN::Neg2LogLikelihood_FCN( const Experiment* data, const Prediction* pdf, const double alpha ) :
  _data( data ), _pdf( pdf ) {

  _pars.Add( "alpha", alpha, 2.e-4, 0., 16. );

}

Neg2LogLikelihood_FCN::~Neg2LogLikelihood_FCN() {
}

double Neg2LogLikelihood_FCN::operator()() const {

  return ( *this )( _pars.Params() );

}

double Neg2LogLikelihood_FCN::operator()( const std::vector< double >& par ) const {

  const Prediction& pdf = *_pdf;
  double result = 0.;
  foreach( const double& mjj, _data->mjjs() )
        {
          const MjjExperiment& mjjData = ( *_data )[mjj];
          for( int bin = 0; bin < mjjData.chi().size(); ++bin ) {
            double x = mjjData.chi( bin );
            int n = mjjData.n( bin );
            double logProb = pdf( mjj, x, n, par );
            result += -2 * logProb;
          }
        }

  return result;
}

double Neg2LogLikelihood_FCN::Up() const {
  return 1.;
}

void Neg2LogLikelihood_FCN::data( const Experiment* data ) {
  _data = data;
  _isMinimized = false;
}

void Neg2LogLikelihood_FCN::pdf( const Prediction* pdf ) {
  _pdf = pdf;
  _isMinimized = false;
}

const Experiment* Neg2LogLikelihood_FCN::data() const {
  return _data;
}
const Prediction* Neg2LogLikelihood_FCN::pdf() const {
  return _pdf;
}
bool Neg2LogLikelihood_FCN::isMinimized() const {
  return _isMinimized;
}

double Neg2LogLikelihood_FCN::Minimize() {
  return Minimize( _pars );
}

double Neg2LogLikelihood_FCN::Minimize( MnUserParameters& pars ) {

  //  cout << "Experiment: " << _data->name() << '\n';
  //  cout << "start value L( " << pars.Params().at( 0 ) << " ) = " << ( *this )( pars.Params() ) << '\n';

  bool isFirst = true;
  vector< double > initVal = pars.Params();
  CombinedMinimizer combMin;
  MnStrategy strat( 2 );
  FunctionMinimum minResults( combMin.Minimize( *this, pars, strat, 10000, 0.1 ) );

  _pars = pars = minResults.UserParameters();
  _isMinimized = minResults.IsValid();

  //  cout << "final value L( " << pars.Params().at( 0 ) << " ) = " << ( *this )( pars.Params() ) << '\n';
  if( !_isMinimized ) cout << "invalid fit!\n";
  return pars.Params().at( 0 );

}

vector< double > Neg2LogLikelihood_FCN::pars() const {
  return _pars.Params();
}

void Neg2LogLikelihood_FCN::pars( const vector< double >& p ) {

  for( int i = 0; i < _pars.Params().size() && i < p.size(); ++i ) {
    //    cout << "setting " << _pars.GetName( i ) << " to: " << p.at( i ) << '\n';
    _pars.SetValue( i, p.at( i ) );
  }

}

void Neg2LogLikelihood_FCN::accept( TestStatMonitor& mon ) {
  mon.monitor( *this );
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

Neg2LogLikelihoodRatio::Neg2LogLikelihoodRatio() :
  _data( 0 ) {
}

Neg2LogLikelihoodRatio::Neg2LogLikelihoodRatio( const Experiment* data, const Prediction* pdf, const double& alpha ) :
  _data( data ), _pdf( pdf ), _numerator( data, pdf, alpha ), _denominator( data, pdf, alpha ) {

}

Neg2LogLikelihoodRatio::~Neg2LogLikelihoodRatio() {
  delete _pdf;
}

double Neg2LogLikelihoodRatio::Up() const {
  return 1.;
}

void Neg2LogLikelihoodRatio::data( const Experiment* data ) {
  _data = data;
  _numerator.data( data );
  _denominator.data( data );
}

void Neg2LogLikelihoodRatio::pdf( const Prediction* pdf ) {
  _pdf = pdf;
  _numerator.pdf( pdf );
  _denominator.pdf( pdf );
}

const Experiment* Neg2LogLikelihoodRatio::data() const {
  return _data;
}

const Prediction* Neg2LogLikelihoodRatio::pdf() const {
  return _pdf;
}

void Neg2LogLikelihoodRatio::accept( TestStatMonitor& mon ) {
  mon.monitor( *this );
}

Neg2LogLikelihood_FCN Neg2LogLikelihoodRatio::numerator() const {
  return _numerator;
}

Neg2LogLikelihood_FCN Neg2LogLikelihoodRatio::denominator() const {
  return _denominator;
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
Neg2LogMaximumLikelihoodRatio::Neg2LogMaximumLikelihoodRatio() :
  Neg2LogLikelihoodRatio(), _inCaseShitMonitor( 0 ) {
}

Neg2LogMaximumLikelihoodRatio::Neg2LogMaximumLikelihoodRatio( const Experiment* data, const Prediction* pdf,
                                                              const double& alpha ) :
  Neg2LogLikelihoodRatio( data, pdf, alpha ), _inCaseShitMonitor( 0 ) {

  init();

}

Neg2LogMaximumLikelihoodRatio::~Neg2LogMaximumLikelihoodRatio() {
  delete _inCaseShitMonitor;
}

void Neg2LogMaximumLikelihoodRatio::init() {
  _denominator.Minimize();
}

double Neg2LogMaximumLikelihoodRatio::operator()( const std::vector< double >& par ) {

  // would nomlize denominator over nuisance parameters .. but none for now
  if( !_denominator.isMinimized() ) _denominator.Minimize();
  if( !_denominator.isMinimized() || _numerator( par ) - _denominator() < 0. ) {

    // try to re-minimize to new value
    _denominator.pars( par );
    _denominator.Minimize();

    // problem?
    if( _numerator( par ) - _denominator() < -0.01 ) {
      // map out the likelihood
      _inCaseShitMonitor = new TestStatMonitor( par.at( 0 ), "figures/", ".png" );
      _numerator.accept( *_inCaseShitMonitor );
      _denominator.accept( *_inCaseShitMonitor );
      _inCaseShitMonitor->finalize();

      cout << "alpha = " << par.at( 0 ) << "\noptimized = " << _denominator.pars().at( 0 ) << "\n numerator   = "
          << _numerator( par ) << "\n denominator = " << _denominator() << "\n ratio       = " << _numerator( par )
          - _denominator() << '\n';

      // _pdf->plot( *_data ); make an experiment visitor for this
      cout << *_data << endl;

      // throw up
      throw runtime_error(
                           lexical_cast< string > ( __FILE__ ) + " line " + lexical_cast< string > ( __LINE__ )
                               + ": Could not minimize denominator." );
    }
  }

  return _numerator( par ) - _denominator();

}

void Neg2LogMaximumLikelihoodRatio::data( const Experiment* data ) {
  _data = data;
  _numerator.data( data );
  _denominator.data( data );
  init();
}

void Neg2LogMaximumLikelihoodRatio::pdf( const Prediction* pdf ) {
  _pdf = pdf;
  _numerator.pdf( pdf );
  _denominator.pdf( pdf );
  init();
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
Neg2LogSimpleLikelihoodRatio::Neg2LogSimpleLikelihoodRatio() :
  Neg2LogLikelihoodRatio() {

}

Neg2LogSimpleLikelihoodRatio::Neg2LogSimpleLikelihoodRatio( const Experiment* data, const Prediction* pdf,
                                                            const double& alpha ) :
  Neg2LogLikelihoodRatio( data, pdf, alpha ) {

  _denominator.pars( vector< double > ( 1, 0. ) );

}

Neg2LogSimpleLikelihoodRatio::~Neg2LogSimpleLikelihoodRatio() {

}

double Neg2LogSimpleLikelihoodRatio::operator()( const std::vector< double >& par ) {
  return _numerator( par ) - _denominator();
}
