#include "Likelihood.hpp"

#include <cfloat>

#include <TH1.h>
#include <TFitterMinuit.h>

#include "PDF.hpp"

using namespace std;


Likelihood::Likelihood()
  : _data( 0 )
{}


Likelihood::Likelihood( const TH1* data, const PDF* pdf)
  : _data( data )
  , _pdf ( pdf  )
{}


Likelihood::~Likelihood() {
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


//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------


LikelihoodRatio::LikelihoodRatio()
  : _data( 0 ) {
}


LikelihoodRatio::LikelihoodRatio( TFitterMinuit * fitter, const TH1* data, const PDF* pdf )
  : _data( (TH1*)data->Clone("data") )
  , _pdf ( pdf )
  , _numerator( data, pdf )
  , _denominator ( 1. )
  , _denominatorL( data, pdf )
  , _fitter( fitter ) {

  _fitter->SetMinuitFCN( &_denominatorL );
  _fitter->SetMaxIterations(5000);
  _fitter->SetPrintLevel(0);
  _fitter->SetStrategy(2);
  _fitter->CreateMinimizer();
  init();

}


LikelihoodRatio::~LikelihoodRatio() {
  _data->Delete();
}


void
LikelihoodRatio::init() {

  vector<double> vec( 1, 0. );
  _fitter->SetParameter( 0, "alpha", vec[0], 1.e-6, 0., 4.e-6);
  int nError = _fitter->Minimize();
  if (nError != 0) {
    cout << ": Problem with Minimize(): " << nError << endl;
    throw( "Problem with Minimization" );
  }

  vec.at(0) = _fitter->GetParameter(0);
  _denominator = _denominatorL( vec );

  // cout << " optimized alpha = " << vec.at(0) << endl;
  // cout << " optimized -2lnL = " << _denominator << endl;

}

double
LikelihoodRatio::operator() ( const std::vector<double>& par ) const {

  // would nomlize denominator over nuisance parameters .. but none for now
  if ( _numerator( par ) < _denominator ) {
    cout << "given alpha        = " << par.at(0) << '\n'
 	 << "optimized alpha    = " << _fitter->GetParameter(0) << '\n'
	 << "-2log(numerator)   = " << _numerator( par ) << '\n'
	 << "-2log(denominator) = " << _denominator      << '\n'
	 << "numerator / denominator = " << exp(-0.5*_numerator( par ) ) << " / " << exp(-0.5*_denominator ) << endl;
  }
  return _numerator( par ) - _denominator ;

}


double LikelihoodRatio::Up() const { return 1.; }

void LikelihoodRatio::data( const TH1* data ) {
  _data = (TH1*)data->Clone("data");
  _numerator.   data( data );
  _denominatorL.data( data );
  init();
}


void LikelihoodRatio::pdf ( const PDF* pdf )   {
  _pdf  = pdf;
  _numerator.   pdf( pdf );
  _denominatorL.pdf( pdf );
  init();
}

TH1* LikelihoodRatio::data() const { return _data; }
const PDF* LikelihoodRatio::pdf () const { return _pdf; }
