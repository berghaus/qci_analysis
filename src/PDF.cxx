#include "PDF.hpp"

#include <cfloat>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include "PDFMonitor.hpp"

using namespace std;

#define foreach BOOST_FOREACH
using boost::lexical_cast;

PDF::PDF() :
    _nData( 1 ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 4. ) ) {
}

PDF::PDF( TFile* file, const double nData ) :
    _nData( nData ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 4. ) ) {

  if ( !file ) throw domain_error(
      lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ ) + ": No input file for PDF" );

  TDirectory * dir = (TDirectory*) file->Get( "2000-mjj-7000GeV" );
  TList * list = dir->GetListOfKeys();
  TIter nextkey( list );
  TKey * key = 0;

  while ( ( key = (TKey*) nextkey() ) ) {
    string name = key->GetName();

    TObject * obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( "TGraph" ) ) {
      double chi = lexical_cast< double >( name.substr( 7 ) );
      TGraphErrors * graph = (TGraphErrors*) obj;
      _eventCounts[chi] = graph;
      TFitResultPtr fitResult = graph->Fit( _pdfFit, "Q" );

      _pdfFitParams[chi].push_back( _pdfFit->GetParameter( 0 ) );
      _pdfFitParams[chi].push_back( _pdfFit->GetParameter( 1 ) );
      _pdfFitParams[chi].push_back( _pdfFit->GetParameter( 2 ) );
    }
  }

}

PDF::PDF( const map< double, vector< double > >& pdfFitParams, double nData ) :
    _nData( nData ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 4. ) ),
    _pdfFitParams( pdfFitParams ) {
}

PDF::PDF( const PDF& orig ) {

  _pdfFit = orig._pdfFit;
  _nData = orig._nData;

  typedef map< double, TGraphErrors* > chiGraphMap_t;
  foreach( chiGraphMap_t::value_type ec, orig._eventCounts )
    _eventCounts[ec.first] = (TGraphErrors*) ec.second->Clone();

  typedef map< double, vector< double > > chiFitParams_t;
  foreach( chiFitParams_t::value_type params, _pdfFitParams )
    _pdfFitParams[params.first] = params.second;

}

PDF::~PDF() {
  _pdfFit->Delete();
}

void PDF::init() {

}

double PDF::operator()( const double& chi, const int& data, const vector< double >& par ) const {

  if ( par.at( 0 ) != par.at( 0 ) ) return 0.;

  if ( par.size() != 1 ) throw( domain_error( "PDF needs at least the compositeness parameter in passed vector" ) );
  double nMC = interpolate( chi, par.at( 0 ) ); // predicted events at x for parameters

  //cout << "n( alpha = " << par.at(0) << " | Data = " << data << ", Chi = " << chi << " ) : " << nMC << '\n';

  double result = 0;
  if ( nMC <= 0 ) throw logic_error(
      lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ )
      + ": Predicted zero or less events" );
  // return log of poisson probability
  if ( data < 0. ) result = 0.;
  else if ( data == 0.0 ) result = -nMC;
  else result = data * log( nMC ) - nMC - TMath::LnGamma( data + 1. ); // == log(TMath::Poisson( i, p ))

  return result;

}

double PDF::operator()( const double& chi, const double& alpha ) const {
  return interpolate( chi, alpha );
}

double PDF::interpolate( const double& chi, const double& alpha ) const {

  typedef map< double, vector< double > > chiFitParams_t;
  // sum over all chi bins
  double sumOverChi = 0.;
  foreach( chiFitParams_t::value_type ec, _pdfFitParams )
  {
    _pdfFit->SetParameter( 0, ec.second[0] );
    _pdfFit->SetParameter( 1, ec.second[1] );
    _pdfFit->SetParameter( 2, ec.second[2] );
    sumOverChi += _pdfFit->Eval( alpha );
  }
  // find chi value of interest
  double lastDistance = DBL_MAX;
  foreach( chiFitParams_t::value_type ec, _pdfFitParams )
  {
    double distance = fabs( ec.first - chi );
    if ( lastDistance < distance ) break;
    lastDistance = distance;
    _pdfFit->SetParameter( 0, ec.second[0] );
    _pdfFit->SetParameter( 1, ec.second[1] );
    _pdfFit->SetParameter( 2, ec.second[2] );
  }
  return _nData * _pdfFit->Eval( alpha ) / sumOverChi;

}

void PDF::accept( PDFMonitor& mon ) {
  mon.monitor( *this );
}

map< double, TGraphErrors* > PDF::eventCounts() const {
  return _eventCounts;
}

map< double, vector< double > > PDF::pdfFitParams() const {
  return _pdfFitParams;
}

TF1* PDF::pdfFit( const std::string& name ) const {
  return (TF1*) _pdfFit->Clone( name.c_str() );
}
