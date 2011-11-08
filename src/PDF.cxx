#include "PDF.hpp"

#include <iostream>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2.h>
#include <TGraph2D.h>

#include "PDFMonitor.hpp"

using namespace std;
using boost::lexical_cast;

PDF::PDF() :
    _hist( 0 ),
    _graph( 0 ) {
}

PDF::PDF( TH2* hist ) :
    _hist( (TH2*) hist->Clone( "PDF" ) ),
    _graph( new TGraph2D( hist ) ) {
}

PDF::PDF( TGraph2D* graph ) :
    _hist( graph->GetHistogram() ),
    _graph( (TGraph2D*)graph->Clone( "PDFGraph" ) ) {
}

PDF::PDF( const PDF& orig ) {

  hist( orig.hist() );
  _graph = orig._graph;
}

PDF::~PDF() {
  _hist->Delete();
  _graph->Delete();
}

double PDF::operator()( const double& x, const int& i, const vector< double >& par ) const {

  if ( par.at( 0 ) != par.at( 0 ) ) return 0.;

  if ( par.size() != 1 ) throw( domain_error( "PDF needs at least the compositeness parameter in passed vector" ) );
  double p = interpolate( x, par.at( 0 ) ); // predicted events at x for parameters

  double result = 0;
  if ( p <= 0 ) throw logic_error(  lexical_cast<string>(__FILE__) +" "+ lexical_cast<string>(__LINE__) +": Predicted zero or less events");
  // return log of poisson probability
  if ( i < 0. ) result = 0.;
  else if ( i == 0.0 ) result = -p;
  else result = i * log( p ) - p - TMath::LnGamma( i + 1. ); // == log(TMath::Poisson( i, p ))

  return result;

}

double PDF::operator()( const double& chi, const double& alpha ) const {
  unsigned int pdfBin = _hist->FindBin( chi, alpha );
  return _hist->GetBinContent( pdfBin );
}

double PDF::interpolate( const double& x, const double& y ) const {
  return _graph->Interpolate( x, y );
}

TH2* PDF::hist() const {
  return _hist;
}

void PDF::hist( const TH2* hist ) {
  _hist = (TH2*) hist->Clone( "PDF" );
}

void PDF::accept( PDFMonitor& mon ) {
  mon.monitor( *this );
}
