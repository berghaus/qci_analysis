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

Likelihood::Likelihood()
  : _data( 0 ) {
  _pars.Add("alpha",0.,2.0e-08,0.,4.e-6);
}


Likelihood::Likelihood( const TH1* data, const PDF* pdf)
  : _data( data )
  , _pdf ( pdf  ) {
  _pars.Add("alpha",0.,2.0e-08,0.,4.e-6);
}


Likelihood::~Likelihood() {
}


double
Likelihood::operator() ( ) const {

  return (*this)( _pars.Params() );

}


double
Likelihood::operator() ( const std::vector<double>& par ) const
{
  const PDF& pdf = *_pdf;

  double result = 0.;
  for ( int bin = 1; bin <= _data->GetNbinsX(); ++bin ) {
    double x = _data->GetBinCenter ( bin );
    int    n = _data->GetBinContent( bin );
    double logProb = pdf( x, n, par );
    result += -2*logProb;
  }

  return result;
}

double
Likelihood::Up() const
{
  return 1.;
}


void Likelihood::data( const TH1* data ) { _data = data; }
void Likelihood::pdf ( const PDF* pdf )  { _pdf  = pdf; }

const TH1* Likelihood::data() const { return _data; }
const PDF* Likelihood::pdf () const { return _pdf; }


double
Likelihood::Minimize(){
  return Minimize( _pars );
}

double
Likelihood::Minimize( MnUserParameters& pars ){

  CombinedMinimizer combMin;
  MnStrategy        strat(2);
  FunctionMinimum   minResults( combMin.Minimize( *this, pars, strat, 10000, 0.1 ) );

  _isMinimized = minResults.IsValid();
  pars = minResults.UserParameters();
  return pars.Params().at(0);
  
}


vector<double> Likelihood::pars() { return _pars.Params(); }


void
Likelihood::accept( TestStatMonitor& monitor ){
}


//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------


LikelihoodRatio::LikelihoodRatio()
  : _data( 0 ) {
}


LikelihoodRatio::LikelihoodRatio( const TH1* data, const PDF* pdf )
  : _data( (TH1*)data->Clone("data") )
  , _pdf ( pdf )
  , _numerator( data, pdf )
  , _denominator( data, pdf ) {

  init();

}


LikelihoodRatio::~LikelihoodRatio() {
  _data->Delete();
}


void
LikelihoodRatio::init() {
  _denominator.Minimize();
}

double
LikelihoodRatio::operator() ( const std::vector<double>& par ) {

  // would nomlize denominator over nuisance parameters .. but none for now
  if ( _numerator( par ) < _denominator() ) {
    cout << "-2lnL("<< par.at(0) << ") = " << _numerator( par ) << '\n'
	 << "-2lnL("<< _denominator.pars().at(0) << ") = " << _denominator() << '\n';
    throw( logic_error("minimized likelihood not at minimum") );
    
  }
  return _numerator( par ) - _denominator() ;

}


double LikelihoodRatio::Up() const { return 1.; }

void LikelihoodRatio::data( const TH1* data ) {
  _data = (TH1*)data->Clone("data");
  _numerator.  data( data );
  _denominator.data( data );
  init();
}


void LikelihoodRatio::pdf ( const PDF* pdf )   {
  _pdf  = pdf;
  _numerator.  pdf( pdf );
  _denominator.pdf( pdf );
  init();
}

TH1* LikelihoodRatio::data() const { return _data; }
const PDF* LikelihoodRatio::pdf () const { return _pdf; }

void
LikelihoodRatio::accept( TestStatMonitor& monitor ){
}
