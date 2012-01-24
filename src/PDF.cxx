#include "PDF.hpp"

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/list_of.hpp>

#include <TClass.h>
#include <TDirectoryFile.h>
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
using namespace boost::assign;

PDF::PDF() :
    _nData( 1 ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 16. ) ) {
}

PDF::PDF( TDirectoryFile* dir, const double nData ) :
    _nData( nData ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 16. ) ) {

  if ( !dir ) throw domain_error(
      lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ ) + ": No input file for PDF" );

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
      TFitResultPtr fitResult = graph->Fit( _pdfFit, "QS" );

      _pdfFitParams[chi].push_back( _pdfFit->GetParameter( 0 ) );
      _pdfFitParams[chi].push_back( _pdfFit->GetParameter( 1 ) );
      _pdfFitParams[chi].push_back( _pdfFit->GetParameter( 2 ) );

      _covarianceMaticies.insert( make_pair( chi, fitResult->GetCovarianceMatrix() ) );

      cout << "insterted at chi = " << chi << endl;

    }
  }

}

PDF::PDF( const map< double, vector< double > >& pdfFitParams, double nData ) :
    _nData( nData ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 4. ) ),
    _pdfFitParams( pdfFitParams ) {
}

PDF::PDF( const PDF& orig ) :
    _pdfFit( static_cast< TF1* >( orig._pdfFit->Clone( "PDFFit" ) ) ),
    _nData( orig._nData ),
    _pdfFitParams( orig.pdfFitParams() ),
    _covarianceMaticies( orig._covarianceMaticies ) {

}

PDF::~PDF() {
  _pdfFit->Delete();
}

void PDF::init() {

}

double PDF::operator()( const double& chi, const int& data, const vector< double >& par ) const {

  if ( par.at( 0 ) != par.at( 0 ) ) return 0.;

  if ( par.size() != 1 ) throw( domain_error( "PDF needs at least the compositeness parameter in passed vector" ) );

  // for looking things up need chi precision rounded to one sig fig
  double x = int( chi * 100 ) % 10 > 4 ?
      ceil( chi * 10 ) / 10. :
      floor( chi * 10 ) / 10.;

  double nMC = interpolate( x, par.at( 0 ) ); // predicted events at x for parameters

  //cout << "n( alpha = " << par.at(0) << " | Data = " << data << ", Chi = " << chi << " ) : " << nMC << " pm " << eMC << '\n';

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
  // for looking things up need chi precision rounded to one sig fig
  double x = int( chi * 100 ) % 10 > 4 ?
      ceil( chi * 10 ) / 10. :
      floor( chi * 10 ) / 10.;

  double nMC = interpolate( x, alpha );
  double eMC = error( x, alpha ); // error on predicted value

  // predicted number of events modified by statistical uncertainty
  // nMC = _random.Gaus( nMC, eMC ); // must be positive

  return nMC;
}

double PDF::interpolate( const double& chi, const double& alpha ) const {

  typedef map< double, vector< double > > chiFitParams_t;
  // sum over all chi bins
  double chiSum = sumOverChi( alpha );

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
  return _nData * _pdfFit->Eval( alpha ) / chiSum;

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

double PDF::sumOverChi( const double& alpha ) const {
  typedef map< double, vector< double > > chiFitParams_t;
  double result = 0.;
  foreach( const chiFitParams_t::value_type ec, _pdfFitParams )
  {
    _pdfFit->SetParameter( 0, ec.second[0] );
    _pdfFit->SetParameter( 1, ec.second[1] );
    _pdfFit->SetParameter( 2, ec.second[2] );
    result += _pdfFit->Eval( alpha );
  }
  return result;
}

double PDF::error( const double& chi, const double& alpha ) const {

  if ( _covarianceMaticies.find( chi ) == _covarianceMaticies.end() ) return 0;

  const TMatrixTSym< double > & mat = _covarianceMaticies.find( chi )->second;

  // I know the gradients of my fit function are:
  vector<double> grad = list_of(0.)( alpha )( pow( alpha, 0.5 ) );
  vector<double> sum(3,0.);
  double error = 0.;

  // first sum
  for( int i = 0; i < mat.GetNrows(); ++i ) {
    sum[i] = 0;
    for( int j = 0; j < mat.GetNcols(); ++j ) {
      sum[i] += mat( i, j ) * grad[j];
    }
  }

  // second sum
  for( int i = 0; i < mat.GetNcols(); ++i ) {
    error += sum[i] * grad[i];
  }

  return error;

}
