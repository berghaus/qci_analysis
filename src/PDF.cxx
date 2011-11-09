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

#include "PDFMonitor.hpp"

using namespace std;

#define foreach BOOST_FOREACH
using boost::lexical_cast;

PDF::PDF() :
    _file( 0 ) {
}

PDF::PDF( TFile* file ) :
    _file( file ) {
  init();
}

PDF::PDF( const PDF& orig ) {
  typedef map< double, TGraphErrors* > chiGraphMap_t;
  foreach( chiGraphMap_t::value_type ec, orig._eventCounts )
    _eventCounts[ec.first] = (TGraphErrors*) ec.second->Clone();
}

PDF::~PDF() {
  _file->Close();
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
      _eventCounts[chi] = (TGraphErrors*) obj;
    }
  }
}

double PDF::operator()( const double& x, const int& i, const vector< double >& par ) const {

  if ( par.at( 0 ) != par.at( 0 ) ) return 0.;

  if ( par.size() != 1 ) throw( domain_error( "PDF needs at least the compositeness parameter in passed vector" ) );
  double p = interpolate( x, par.at( 0 ) ); // predicted events at x for parameters

  double result = 0;
  if ( p <= 0 ) throw logic_error(
      lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ )
      + ": Predicted zero or less events" );
  // return log of poisson probability
  if ( i < 0. ) result = 0.;
  else if ( i == 0.0 ) result = -p;
  else result = i * log( p ) - p - TMath::LnGamma( i + 1. ); // == log(TMath::Poisson( i, p ))

  return result;

}

double PDF::operator()( const double& chi, const double& alpha ) const {
  return interpolate( chi, alpha );
}

double PDF::interpolate( const double& chi, const double& alpha ) const {

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

void PDF::accept( PDFMonitor& mon ) {
  mon.monitor( *this );
}

map<double,TGraphErrors*> PDF::eventCounts() const {
  typedef map< double, TGraphErrors* > chiGraphMap_t;
  chiGraphMap_t result;
  foreach( chiGraphMap_t::value_type ec, _eventCounts )
    result[ec.first] = (TGraphErrors*) ec.second->Clone();
  return result;
}
