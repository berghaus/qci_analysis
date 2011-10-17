#include "PDF.hpp"

#include <iostream>
#include <stdexcept>

#include <TMath.h>
#include <TH2.h>

using namespace std;

PDF::PDF()
  : _hist( 0 ) {
}



PDF::PDF( const TH2* hist ) 
  : _hist ( (TH2*)hist->Clone("PDF") ) {
}


PDF::PDF( const PDF& orig ) {

  hist( orig.hist() );

}


PDF::~PDF() {
  cerr << "PDF destructor" << endl;
  _hist->Delete();
}


double
PDF::operator() (  const double& x, const int& i, const vector<double>& par ) const{

  if ( par.size() != 1 ) throw( domain_error("PDF needs at least the compositeness parameter in passed vector") );
  double p = _hist->Interpolate( x, par.at(0) ); // predicted events at x for parameters

  return TMath::Poisson( i, p );   // probability of observing i having predicted p

}


TH2* PDF::hist() const { return _hist; }


void PDF::hist( const TH2* hist ) { _hist = (TH2*)hist->Clone("PDF"); }
