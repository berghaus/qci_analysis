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

#include "PDFMonitor.hpp"

using namespace std;

#define foreach BOOST_FOREACH
using boost::lexical_cast;

PDF::PDF() :
    _file( 0 ),
    _nData( 1 ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x", 0., 4. ) ),
    _normalizedPdfFit( new TF1( "PDFFit", "([0]+[1]*x)/([2]+[3]*x)*[4]", 0., 4. ) ),
    _useFit( 0 ) {
}

PDF::PDF( TFile* file, const double nData ) :
    _file( file ),
    _nData( nData ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x", 0., 4. ) ),
    _normalizedPdfFit( new TF1( "PDFFit", "([0]+[1]*x)/([2]+[3]*x)*[4]", 0., 4. ) ),
    _useFit( 0 ) {
  init();
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

  _useFit = orig._useFit;

}

PDF::~PDF() {
  _pdfFit->Delete();
  _normalizedPdfFit->Delete();
}

void PDF::init() {
  if ( !_file ) throw domain_error(
      lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ ) + ": No input file for PDF" );
  TDirectory * dir = (TDirectory*) _file->Get( "2000-mjj-7000GeV" );
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
      for( int i = 0; i < graph->GetN(); ++i ) {
        cout << "(n, alpha, n(alpha) ) = ( " << i << ", " << graph->GetX()[i] << ", " << graph->GetY()[i] << " )\n";
        if ( graph->GetX()[i] == 0. ) _pdfFit->FixParameter(0, graph->GetY()[i] );
      }
      graph->Fit( _pdfFit );
      _pdfFitParams[chi].push_back( _pdfFit->GetParameter( 0 ) );
      _pdfFitParams[chi].push_back( _pdfFit->GetParameter( 1 ) );
      if ( chi == 0. ) {
        _normalizedPdfFit->SetParameter( 2, _pdfFit->GetParameter( 0 ) );
        _normalizedPdfFit->SetParameter( 3, _pdfFit->GetParameter( 1 ) );
        _normalizedPdfFit->SetParameter( 4, _nData );
      }
    }
  }
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

  if ( _useFit ) {
    typedef map< double, vector< double > > chiFitParams_t;
    double lastDistance = DBL_MAX;
    foreach( chiFitParams_t::value_type ec, _pdfFitParams )
    {
      double distance = fabs( ec.first - chi );
      if ( lastDistance < distance ) break;
      lastDistance = distance;
      _pdfFit->SetParameter( 0, ec.second[0] );
      _pdfFit->SetParameter( 1, ec.second[1] );

      _normalizedPdfFit->SetParameter( 0, ec.second[0] );
      _normalizedPdfFit->SetParameter( 1, ec.second[1] );
    }
    return _normalizedPdfFit->Eval( alpha );

  } else {

    // determine closest TGraphErrors
    typedef map< double, TGraphErrors* > chiGraphMap_t;
    const chiGraphMap_t::value_type * relevantEntry = 0;
    double lastDistance = DBL_MAX;
    foreach( const chiGraphMap_t::value_type& ec, _eventCounts )
    {
      double distance = fabs( ec.first - chi );
      if ( lastDistance < distance ) break;
      lastDistance = distance;
      relevantEntry = &ec;
    }
    return relevantEntry->second->Eval( alpha );

  }

  return 0.;

}

void PDF::accept( PDFMonitor& mon ) {
  mon.monitor( *this );
}

map< double, TGraphErrors* > PDF::eventCounts() const {
  typedef map< double, TGraphErrors* > chiGraphMap_t;
  chiGraphMap_t result;
  foreach( chiGraphMap_t::value_type ec, _eventCounts )
    result[ec.first] = (TGraphErrors*) ec.second->Clone();
  return result;
}

map< double, vector< double > > PDF::pdfFitParams() const {
  return _pdfFitParams;
}

void PDF::useFit( const bool useFit ) {
  _useFit = useFit;
}

TF1* PDF::pdfFit( const std::string& name ) const {
  return (TF1*) _pdfFit->Clone( name.c_str() );
}
