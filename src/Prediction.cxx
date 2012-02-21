#include "Prediction.hpp"

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

#include "PredictionMonitor.hpp"

using namespace std;

#define foreach BOOST_FOREACH
using boost::lexical_cast;
using namespace boost::assign;

//_____________________________________________________________________________________________________________________
Prediction::Prediction() :
    _nData( 1 ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 16. ) ) {
}

//_____________________________________________________________________________________________________________________
Prediction::Prediction( TDirectoryFile* dir, const double nData ) :
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

//_____________________________________________________________________________________________________________________
Prediction::Prediction( const map< double, vector< double > >& pdfFitParams, double nData ) :
    _nData( nData ),
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 4. ) ),
    _pdfFitParams( pdfFitParams ) {
}

//_____________________________________________________________________________________________________________________
Prediction::Prediction( const Prediction& orig ) :
    _pdfFit( static_cast< TF1* >( orig._pdfFit->Clone( "PDFFit" ) ) ),
    _nData( orig._nData ),
    _pdfFitParams( orig.pdfFitParams() ),
    _covarianceMaticies( orig._covarianceMaticies ) {

}

//_____________________________________________________________________________________________________________________
Prediction::~Prediction() {
  _pdfFit->Delete();
}

//_____________________________________________________________________________________________________________________
void Prediction::init() {
}

//_____________________________________________________________________________________________________________________
double Prediction::operator()( const double& chi, const int& data, const vector< double >& par ) const {

  if ( par.at( 0 ) != par.at( 0 ) ) return 0.;

  if ( par.size() != 1 ) throw( domain_error( "PDF needs at least the compositeness parameter in passed vector" ) );

  // for looking things up need chi precision rounded to one sig fig
  double x = labelChi( chi );

  // predicted events at chi for parameters
  double nMC = interpolate( x, par.at( 0 ) );

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

//_____________________________________________________________________________________________________________________
double Prediction::operator()( const double& chi, const double& alpha ) const {

  // for looking things up need chi precision rounded to one sig fig
  double x = labelChi( chi );
  double mjj = 2001.;
  double lambda = pow( alpha, -0.25 );

  double nMC = interpolate( x, alpha );

//  cout << "\n--------------------------------\n";
//  cout << "starting nMC = " << nMC << "\n";
  foreach( const Effect* eff, _effects ) nMC = eff->apply( nMC, lambda, x, mjj );
//  cout << "final nMC = " << nMC << "\n";
//  cout << "outside error = " << eMC << endl;
//  cout << "--------------------------------\n";

  // predicted number of events modified by statistical uncertainty
  //nMC = _random.Gaus( nMC, eMC ); // must be positive

  return nMC;
}

//_____________________________________________________________________________________________________________________
double Prediction::interpolate( const double& chi, const double& alpha ) const {

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

//_____________________________________________________________________________________________________________________
void Prediction::accept( PredictionMonitor& mon ) {
  mon.monitor( *this );
}

//_____________________________________________________________________________________________________________________
map< double, TGraphErrors* > Prediction::eventCounts() const {
  return _eventCounts;
}

//_____________________________________________________________________________________________________________________
map< double, vector< double > > Prediction::pdfFitParams() const {
  return _pdfFitParams;
}

//_____________________________________________________________________________________________________________________
TF1* Prediction::pdfFit( const std::string& name ) const {
  return (TF1*) _pdfFit->Clone( name.c_str() );
}

//_____________________________________________________________________________________________________________________
double Prediction::sumOverChi( const double& alpha ) const {
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

//_____________________________________________________________________________________________________________________
void Prediction::nData( const int& nData ) {
  _nData = nData;
}

//_____________________________________________________________________________________________________________________
int Prediction::nData() const {
  return _nData;
}

//_____________________________________________________________________________________________________________________
double Prediction::labelChi( const double& chi ) const {
  return int( chi * 100 ) % 10 > 4 ?
      ceil( chi * 10 ) / 10. :
      floor( chi * 10 ) / 10.;
}

//_____________________________________________________________________________________________________________________
map< double, TMatrixTSym< double > > Prediction::covarianceMaticies() const {
  return _covarianceMaticies;
}

//_____________________________________________________________________________________________________________________
void Prediction::effects( vector< Effect* > effects ) {
  _effects = effects;
}

//_____________________________________________________________________________________________________________________
void Prediction::addEffect( Effect* effect ) {
  _effects.push_back( effect );
}

//_____________________________________________________________________________________________________________________
void Prediction::newPE() const {
  foreach( Effect* e, _effects ) e->newPE();
}
