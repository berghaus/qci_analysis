#include "Prediction.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <set>
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
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 16. ) ) {
}

//_____________________________________________________________________________________________________________________
Prediction::Prediction( map< double, TDirectoryFile* > dirs, const Experiment& exp ) :
    _pdfFit( new TF1( "PDFFit", "[0]+[1]*x+[2]*sqrt(x)", 0., 16. ) ) {

  typedef map< double, TDirectoryFile* > dirMap_t;
  foreach( const dirMap_t::value_type mjjDir, dirs )
  {
    double mjj = mjjDir.first;
    TDirectoryFile* dir = mjjDir.second;
    if ( !dir ) throw domain_error(
        lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ ) + ": No input file for PDF" );

    _predictions.insert( make_pair( mjj, MjjPrediction( dir, exp[mjj].integral(), _pdfFit ) ) );
    _mjjs.insert( mjj );

  }

}

//_____________________________________________________________________________________________________________________
Prediction::Prediction( const Prediction& orig ) :
    _pdfFit( static_cast< TF1* >( orig._pdfFit->Clone( "PDFFit" ) ) ),
    _predictions( orig._predictions ),
    _effects( orig._effects ),
    _mjjs( orig._mjjs ) {
}

//_____________________________________________________________________________________________________________________
Prediction::~Prediction() {
  _pdfFit->Delete();
}

//_____________________________________________________________________________________________________________________
void Prediction::init() {
}

//_____________________________________________________________________________________________________________________
double Prediction::operator()( const double& mjj, const double& chi, const int& data,
                               const vector< double >& par ) const {

  if ( par.at( 0 ) != par.at( 0 ) ) return 0.;

  if ( par.size() != 1 ) throw( domain_error( "PDF needs at least the compositeness parameter in passed vector" ) );

  setMjj( mjj );

  // for looking things up need chi precision rounded to one sig fig
  double x = labelChi( chi );

  // predicted events at chi for parameters
  double nMC = interpolate( x, par.at( 0 ) );

  //cout << "n( alpha = " << par.at(0) << " | Data = " << data << ", Chi = " << chi << " ) : " << nMC << " pm " << eMC << '\n';

  double result = 0;
  if ( nMC <= 0 ) throw logic_error(
      lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ )
      + ": Predicted zero or less events" );

  // return log of Poisson probability
  if ( data < 0. ) result = 0.;
  else if ( data == 0.0 ) result = -nMC;
  else result = data * log( nMC ) - nMC - TMath::LnGamma( data + 1. ); // == log(TMath::Poisson( i, p ))

  return result;

}

//_____________________________________________________________________________________________________________________
double Prediction::operator()( const double& mjj, const double& chi, const double& alpha ) const {

  // for looking things up need chi precision rounded to one sig fig
  setMjj( mjj );
  double x = labelChi( chi );
  double lambda = pow( alpha, -0.25 );

  double nMC = interpolate( x, alpha );

//  cout << "\n--------------------------------\n";
//  cout << "starting nMC = " << nMC << "\n";
  foreach( const Effect* eff, _effects )
    nMC = eff->apply( nMC, lambda, x, mjj );
//  cout << "final nMC = " << nMC << "\n";
//  cout << "outside error = " << eMC << endl;
//  cout << "--------------------------------\n";

  return nMC;
}

//_____________________________________________________________________________________________________________________
double Prediction::interpolate( const double& chi, const double& alpha ) const {

  typedef map< double, vector< double > > chiFitParams_t;
  // sum over all chi bins
  double chiSum = sumOverChi( alpha );

  // find chi value of interest
  double lastDistance = DBL_MAX;

  foreach( chiFitParams_t::value_type ec, _mjjPred->_pdfFitParams )
  {
    double distance = fabs( ec.first - chi );
    if ( lastDistance < distance ) break;
    lastDistance = distance;
    _pdfFit->SetParameter( 0, ec.second[0] );
    _pdfFit->SetParameter( 1, ec.second[1] );
    _pdfFit->SetParameter( 2, ec.second[2] );
  }
  return _mjjPred->_nData * _pdfFit->Eval( alpha ) / chiSum;

}

//_____________________________________________________________________________________________________________________
void Prediction::accept( PredictionMonitor& mon ) {
  mon.monitor( *this );
}

//_____________________________________________________________________________________________________________________
map< double, TGraphErrors* > Prediction::eventCounts( const double& mjj ) const {
  setMjj( mjj );
  return _mjjPred->_eventCounts;
}

//_____________________________________________________________________________________________________________________
map< double, vector< double > > Prediction::pdfFitParams( const double& mjj ) const {
  setMjj( mjj );
  return _mjjPred->_pdfFitParams;
}

//_____________________________________________________________________________________________________________________
TF1* Prediction::pdfFit( const std::string& name ) const {
  return (TF1*) _pdfFit->Clone( name.c_str() );
}

//_____________________________________________________________________________________________________________________
double Prediction::sumOverChi( const double& alpha ) const {
  typedef map< double, vector< double > > chiFitParams_t;
  double result = 0.;
  foreach( const chiFitParams_t::value_type ec, _mjjPred->_pdfFitParams )
  {
    _pdfFit->SetParameter( 0, ec.second[0] );
    _pdfFit->SetParameter( 1, ec.second[1] );
    _pdfFit->SetParameter( 2, ec.second[2] );
    result += _pdfFit->Eval( alpha );
  }
  return result;
}

//_____________________________________________________________________________________________________________________
void Prediction::nData( const Experiment& e ) {

  set<double> eMjjs = e.mjjs();
  vector<double> dummy( max( _mjjs.size(), eMjjs.size() ) );
  vector<double>::iterator it = dummy.begin();
  it = set_intersection( _mjjs.begin(), _mjjs.end(), eMjjs.begin(), eMjjs.end(), it );
  vector<double> overlap;
  overlap.insert( overlap.begin(), dummy.begin(), it );
  if ( overlap.empty() ) {
    cout << "experiment mjj = [ ";
    foreach( double mjj, eMjjs ) cout << mjj << " ";
    cout << "]\n";

    cout << "prediction mjj = [ ";
    foreach( double mjj, _mjjs ) cout << mjj << " ";
    cout << "]\n";

    cout << "intersection mjj = [ ";
    foreach( double mjj, overlap ) cout << mjj << " ";
    cout << "]\n";

    throw domain_error(
      lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ ) + ": No overlap in m_jj " );
  }
  foreach( const double& mjj, overlap ) {
    setMjj( mjj );
    _mjjPred->_nData = e[mjj].integral();
  }

}

//_____________________________________________________________________________________________________________________
int Prediction::nData( const double& mjj ) const {
  setMjj( mjj );
  return _mjjPred->_nData;
}

//_____________________________________________________________________________________________________________________
double Prediction::labelChi( const double& chi ) const {
  return int( chi * 100 ) % 10 > 4 ?
      ceil( chi * 10 ) / 10. :
      floor( chi * 10 ) / 10.;
}

//_____________________________________________________________________________________________________________________
map< double, TMatrixTSym< double > > Prediction::covarianceMaticies( const double& mjj ) const {
  setMjj( mjj );
  return _mjjPred->_covarianceMaticies;
}

//_____________________________________________________________________________________________________________________
map< double, map< double, TMatrixTSym< double > > > Prediction::covarianceMaticies() const {
  map< double, map< double, TMatrixTSym< double > > > result;
  foreach( const double& mjj, _mjjs )
    result.insert( make_pair( mjj, covarianceMaticies( mjj ) ) );
  return result;
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
  foreach( Effect* e, _effects )
    e->newPE();
}

//_____________________________________________________________________________________________________________________
void Prediction::setMjj( const double& mjj ) const {
  _mjj = mjj;
  std::map< double, MjjPrediction >::const_iterator itr = _predictions.lower_bound( _mjj );
  if ( itr == _predictions.end() ) throw domain_error(
      lexical_cast< string >( __FILE__ ) + " " + lexical_cast< string >( __LINE__ ) + ": No Prediction for mjj = "
      + lexical_cast< string >( _mjj ) );
  _mjjPred = const_cast< MjjPrediction* >( &( itr->second ) );
}


//_____________________________________________________________________________________________________________________
Prediction::MjjPrediction::MjjPrediction( TDirectoryFile* dir, const double nData, TF1* pdfFit ) :
    _nData( nData ),
    _pdfFit( pdfFit ) {

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
