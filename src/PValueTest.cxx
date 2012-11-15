#include "PValueTest.hpp"
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/checked_delete.hpp>

#include <TCanvas.h>
#include <TFitterMinuit.h>
#include <TH1.h>
#include <TArrow.h>
#include <TText.h>

#include "Prediction.hpp"
#include "Likelihood.hpp"

#define foreach BOOST_FOREACH
using namespace std;
using boost::format;


template< typename Likelihood >
PValueTest<Likelihood>::PValueTest( const double alpha, const vector< Likelihood* >& lambdas ) :
    _alpha( alpha ),
    _lambdas( lambdas ),
    _dataLLR( 0 ) {
  init();
}


template< typename Likelihood >
void PValueTest<Likelihood>::init() {

  vector< double > par( 1, _alpha );
  _testStats.reserve( _lambdas.size() );
  foreach( Likelihood* lambda, _lambdas )
  {
    Likelihood& l = *lambda;
    _testStats.push_back( l( par ) );
  }
  sort( _testStats.begin(), _testStats.end() );

}

template< typename Likelihood >
PValueTest<Likelihood>::PValueTest() :
    _alpha( 0 ),
    _dataLLR( 0 ) {
}

template< typename Likelihood >
PValueTest<Likelihood>::~PValueTest() {
}

template< typename Likelihood >
double PValueTest<Likelihood>::operator()( Likelihood& lambda ) {

  vector< double > par( 1, _alpha );

  return ( *this )( lambda( par ) );

}

template< typename Likelihood >
double PValueTest<Likelihood>::operator()( const double& llr ) {

  _dataLLR = llr;
  vector< double >::iterator itr = lower_bound( _testStats.begin(), _testStats.end(), _dataLLR );
  int nOver = distance( itr, _testStats.end() );

  return double( nOver ) / double( _testStats.size() );

}

template< typename Likelihood >
void PValueTest<Likelihood>::finalize( const std::string& dir ) {

  double histMax = _testStats[_testStats.size() / 2] > _dataLLR ?
      2 * _testStats[_testStats.size() / 2] :
      2 * _dataLLR;
  if ( histMax < 1. ) histMax = 1.;
  int nBins = _testStats.size() / 100;
  //double off = 1.5 * histMax / double( nBins );
  double off = 0;

  string hName = str( format( "Likelihood_FCN-scale%2.2e" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(), "", nBins, -off, histMax - off );
  string xTitle = _alpha == 0 ? "-2Likelihood( #Lambda = #infty TeV )"
                              : str( format( "Likelihood( #Lambda = %2.2f TeV )" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );

  foreach( const double& x, _testStats )
    _minus2LnLikelihoodDistribution->Fill( x );


  double middle = (_minus2LnLikelihoodDistribution->GetMaximum() - _minus2LnLikelihoodDistribution->GetMinimum())/2;
  TArrow* dataArrow = new TArrow( _dataLLR / 2, middle, _dataLLR*.99, middle*0.02  );

  TText * dataText = new TText( _dataLLR / 3, middle*1.05,  "Data" );
  dataArrow->SetLineWidth( 2 );

  TCanvas* pvc = new TCanvas( ( hName + "Canvas" ).c_str(), "", 500, 500 );
  pvc->cd();
  _minus2LnLikelihoodDistribution->Draw();
  int dataBin = _minus2LnLikelihoodDistribution->FindBin( _dataLLR );
  int lastBin = _minus2LnLikelihoodDistribution->GetNbinsX();
  TH1 * ldClone = (TH1*) _minus2LnLikelihoodDistribution->Clone();
  ldClone->SetFillColor( kRed );
  ldClone->SetFillStyle( 1001 );
  ldClone->GetXaxis()->SetRange( dataBin, lastBin );
  ldClone->Draw("SAME][");
  dataArrow->Draw();
  dataText->Draw();
  string cName = _minus2LnLikelihoodDistribution->GetName();
  pvc->Print( ( dir + cName + ".pdf" ).c_str() );
}

template<>
void PValueTest<Neg2LogLikelihood_FCN>::finalize( const std::string& dir ) {

  double histMax = _testStats[_testStats.size() / 2] > _dataLLR ?
      2 * _testStats[_testStats.size() / 2] :
      2 * _dataLLR;
  if ( histMax < 1. ) histMax = 1.;
  int nBins = _testStats.size() / 100;
  //double off = 1.5 * histMax / double( nBins );
  double off = 0;

  string hName = str( format( "Likelihood_FCN-scale%2.2e" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(), "", nBins, -off, histMax - off );
  string xTitle = _alpha == 0 ? "-2lnL( n | #Lambda = #infty TeV )"
                              : str( format( "-2lnL( n | #Lambda = %2.2f TeV )" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );

  foreach( const double& x, _testStats )
    _minus2LnLikelihoodDistribution->Fill( x );


  double middle = (_minus2LnLikelihoodDistribution->GetMaximum() - _minus2LnLikelihoodDistribution->GetMinimum())/2;
  TArrow* dataArrow = new TArrow( _dataLLR / 2, middle, _dataLLR*.99, middle*0.02  );

  TText * dataText = new TText( _dataLLR / 3, middle*1.05,  "Data" );
  dataArrow->SetLineWidth( 2 );

  TCanvas* pvc = new TCanvas( ( hName + "Canvas" ).c_str(), "", 500, 500 );
  pvc->cd();
  _minus2LnLikelihoodDistribution->Draw();
  int dataBin = _minus2LnLikelihoodDistribution->FindBin( _dataLLR );
  int lastBin = _minus2LnLikelihoodDistribution->GetNbinsX();
  TH1 * ldClone = (TH1*) _minus2LnLikelihoodDistribution->Clone();
  ldClone->SetFillColor( kRed );
  ldClone->SetFillStyle( 1001 );
  ldClone->GetXaxis()->SetRange( dataBin, lastBin );
  ldClone->Draw("SAME][");
  dataArrow->Draw();
  dataText->Draw();
  string cName = _minus2LnLikelihoodDistribution->GetName();
  pvc->Print( ( dir + cName + ".pdf" ).c_str() );
}

template<>
void PValueTest<Neg2LogMaximumLikelihoodRatio>::finalize( const std::string& dir ) {

//  double histMax = _testStats[_testStats.size() / 2] > _dataLLR ?
//      2 * _testStats[_testStats.size() / 2] :
//      2 * _dataLLR;
  _dataLLR = fabs( _dataLLR );
  double histMax = log10( _testStats.back() );
  int nBins = _testStats.size() / 100;
  double off = -10;

  cout << "dataLLR = " << _dataLLR << endl;
  cout << "test stst front = " << _testStats.front() << " and back = " << _testStats.back() << endl;


  string hName = str( format( "Likelihood_FCN-scale%2.2e" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution = new TH1D( hName.c_str(), "", nBins, off, histMax );
  string xTitle = _alpha == 0 ? "ln( -2lnQ( #Lambda = #infty TeV ) )"
                              : str( format( "log( -2lnQ( #Lambda = %2.2f TeV ) )" ) % pow( _alpha, -0.25 ) );
  _minus2LnLikelihoodDistribution->SetXTitle( xTitle.c_str() );
  _minus2LnLikelihoodDistribution->GetXaxis()->SetNdivisions( 505 );

  foreach( const double& x, _testStats )
    _minus2LnLikelihoodDistribution->Fill( log10(fabs(x)) );


  double middle = (_minus2LnLikelihoodDistribution->GetMaximum() - _minus2LnLikelihoodDistribution->GetMinimum())/10;
  TArrow* dataArrow = new TArrow( (histMax+off) / 2, middle*9.5, log10(_dataLLR) , 1.  );

  TText * dataText = new TText( (histMax+off) / 2,  middle*10,  "Data" );
  dataArrow->SetLineWidth( 2 );

  TCanvas* pvc = new TCanvas( ( hName + "Canvas" ).c_str(), "", 500, 500 );
  pvc->cd();
  pvc->SetLogy();
  _minus2LnLikelihoodDistribution->SetMinimum(0.9);
  _minus2LnLikelihoodDistribution->Draw();
  int dataBin = _minus2LnLikelihoodDistribution->FindBin( log10(_dataLLR) );
  int lastBin = _minus2LnLikelihoodDistribution->GetNbinsX();
  TH1 * ldClone = (TH1*) _minus2LnLikelihoodDistribution->Clone();
  ldClone->SetFillColor( kRed );
  ldClone->SetFillStyle( 1001 );
  ldClone->GetXaxis()->SetRange( dataBin+1, lastBin );
  ldClone->Draw("SAME][");
  dataArrow->Draw();
  dataText->Draw();
  string cName = _minus2LnLikelihoodDistribution->GetName();
  pvc->Print( ( dir + cName + ".pdf" ).c_str() );
}

template< typename Likelihood >
double PValueTest<Likelihood>::alpha() const {
  return _alpha;
}

template< typename Likelihood >
void PValueTest<Likelihood>::alpha( const double& alpha ) {
  _alpha = alpha;
}

template< typename Likelihood >
void PValueTest<Likelihood>::clear() {

  _alpha = -1.;
  _lambdas.clear();
  _dataLLR = 0.;

  //delete _minus2LnLikelihoodDistribution;
  _testStats.clear();

}

template< typename Likelihood >
ostream& operator<<( ostream& out, const PValueTest<Likelihood>& test ) {

  out << test.alpha() << " ";
  vector< double >::const_iterator itr = test._testStats.begin();
  vector< double >::const_iterator end = test._testStats.end();
  for( ; itr != end; ++itr ) {
    out << *itr << " ";
  }

  return out;
}


template< typename Likelihood >
istream& operator>>( istream& in, PValueTest<Likelihood>& test ) {

  test.clear();

  string buffer;
  double x;

  getline( in, buffer );
  if ( !in ) return in;
  stringstream bufferStream( buffer );
  bufferStream >> x;
  test.alpha( x );
  while ( bufferStream >> x ) {
    test._testStats.push_back( x );
  }

  return in;
}


template class PValueTest<Neg2LogLikelihood_FCN>;
template typename std::ostream& operator<< <Neg2LogLikelihood_FCN>( std::ostream&, const PValueTest<Neg2LogLikelihood_FCN>& );
template typename std::istream& operator>> <Neg2LogLikelihood_FCN>( std::istream&, PValueTest<Neg2LogLikelihood_FCN>& );

template class PValueTest<Neg2LogMaximumLikelihoodRatio>;
template typename std::ostream& operator<< <Neg2LogMaximumLikelihoodRatio>( std::ostream&, const PValueTest<Neg2LogMaximumLikelihoodRatio>& );
template typename std::istream& operator>> <Neg2LogMaximumLikelihoodRatio>( std::istream&, PValueTest<Neg2LogMaximumLikelihoodRatio>& );

template class PValueTest<Neg2LogSimpleLikelihoodRatio>;
template typename std::ostream& operator<< <Neg2LogSimpleLikelihoodRatio>( std::ostream&, const PValueTest<Neg2LogSimpleLikelihoodRatio>& );
template typename std::istream& operator>> <Neg2LogSimpleLikelihoodRatio>( std::istream&, PValueTest<Neg2LogSimpleLikelihoodRatio>& );

