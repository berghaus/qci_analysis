#include "PDF.hpp"

#include <iostream>
#include <stdexcept>

#include <TCanvas.h>
#include <TMath.h>
#include <TH2.h>
#include <TGraph2D.h>

#include "PDFMonitor.hpp"

using namespace std;

PDF::PDF() :
    _hist( 0 ),
    _graph( 0 ) {
}

PDF::PDF( TH2* hist ) :
    _hist( (TH2*) hist->Clone( "PDF" ) ),
    _graph( new TGraph2D( hist ) ) {
}

PDF::PDF( TGraph2D* graph ) :
    _hist( _graph->GetHistogram() ),
    _graph( (TGraph2D*)graph->Clone( "PDFGraph" ) ) {
}

PDF::PDF( const PDF& orig ) {

  hist( orig.hist() );
  _graph = new TGraph2D( orig.hist() );
}

PDF::~PDF() {
  _hist->Delete();
  _graph->Delete();
}

double PDF::operator()( const double& x, const int& i, const vector< double >& par ) const {

  if ( par.at( 0 ) != par.at( 0 ) ) return 0.;

  if ( par.size() != 1 ) throw( domain_error( "PDF needs at least the compositeness parameter in passed vector" ) );
  double p = interpolate( x, par.at( 0 ) ); // predicted events at x for parameters

  // return log of poisson probability
  if ( i < 0. ) return 0.;
  else if ( i == 0.0 ) return -p;
  return i * log( p ) - p - TMath::LnGamma( i + 1. ); // == log(TMath::Poisson( i, p ))

}

double PDF::operator()( const double& chi, const double& alpha ) const {
  unsigned int pdfBin = _hist->FindBin( chi, alpha );
  return _hist->GetBinContent( pdfBin );
}

double PDF::interpolate( const double& x, const double& y ) const {
  return _hist->Interpolate( x, y );
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
