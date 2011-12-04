#include "CertaintyLevel.hpp"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH2.h>
#include <TMultiGraph.h>


using namespace std;
using namespace boost::math;

CertaintyLevel::CertaintyLevel() :
  _nBins( 1 ),
  _min( 2. ),
  _max( 8. ) {
}


  // constructor specifying name (CLs, CLs+b, etc) and number of values to expect
CertaintyLevel::CertaintyLevel( const string& name, const vector<double>::size_type& size ) :
  _nBins( size ),
  _min( 2. ),
  _max( 8. ),
  _name( name ) {
  _scales.reserve( size );
  _obs.reserve( size );
  _exp.reserve( size );
  _e02.reserve( size );
  _e98.reserve( size );
  _e16.reserve( size );
  _e84.reserve( size );
  _ex.reserve( size );
}


CertaintyLevel::CertaintyLevel( const string& name, const vector<double>::size_type& size, const double& min, const double& max ) :
  _nBins( size ),
  _min( min ),
  _max( max ),
  _name( name) {
  _scales.reserve( size );
  _obs.reserve( size );
  _exp.reserve( size );
  _e02.reserve( size );
  _e98.reserve( size );
  _e16.reserve( size );
  _e84.reserve( size );
  _ex.reserve( size );
}


void CertaintyLevel::add( const double& scale, const double& obs, const vector<double>& exps ) {

  double sigmaN2 = quantile( exps, cdf( _myGaus, -2. ) );
  double sigmaN1 = quantile( exps, cdf( _myGaus, -1. ) );
  double exp = quantile( exps, cdf( _myGaus, 0. ) );
  double sigmaP1 = quantile( exps, cdf( _myGaus, 1. ) );
  double sigmaP2 = quantile( exps, cdf( _myGaus, 2. ) );

  _scales.push_back( scale );
  _obs.push_back( obs );
  _exp.push_back( exp );
  _e02.push_back( fabs( exp - sigmaN2 ) );
  _e98.push_back( fabs( exp - sigmaP2 ) );
  _e16.push_back( fabs( exp - sigmaN1 ) );
  _e84.push_back( fabs( exp - sigmaP1 ) );
  _ex.push_back( 0.5*(_max-_min)/double(_nBins) );

}


void CertaintyLevel::plot() {

  assert( _scales.size() == _nBins );
  assert( _obs.size() == _nBins );
  assert( _exp.size() == _nBins );
  assert( _e02.size() == _nBins );
  assert( _e98.size() == _nBins );
  assert( _e16.size() == _nBins );
  assert( _e84.size() == _nBins );

  TGraph * observedGraph = new TGraph( _nBins, &_scales[0], &_obs[0] );
  observedGraph->SetLineWidth( 2 );
  observedGraph->SetLineColor( kRed );

  TGraphAsymmErrors * expected1Sigma = new TGraphAsymmErrors( _nBins, &_scales[0], &_exp[0], &_ex[0], &_ex[0], &_e16[0], &_e84[0] );
  expected1Sigma->SetFillStyle( 1001 );
  expected1Sigma->SetFillColor( kGreen );
  expected1Sigma->SetLineStyle( 2 );
  expected1Sigma->SetLineWidth( 2 );

  TGraphAsymmErrors * expected2Sigma = new TGraphAsymmErrors( _nBins, &_scales[0], &_exp[0], &_ex[0], &_ex[0], &_e02[0], &_e98[0] );
  expected2Sigma->SetFillStyle( 1001 );
  expected2Sigma->SetFillColor( kYellow );
  expected2Sigma->SetLineColor( kYellow );
  expected2Sigma->SetLineWidth( 2 );

  TMultiGraph * exclusion = new TMultiGraph;
  exclusion->Add( expected2Sigma, "3" );
  exclusion->Add( expected1Sigma, "3" );
  exclusion->Add( expected1Sigma, "LX" );
  exclusion->Add( observedGraph, "L" );

  TCanvas * fullCanvas = new TCanvas( (_name+"_full").c_str(), "", 800, 800 );
  fullCanvas->cd();
  TH2D * fullHist = new TH2D( (_name+"fullHist").c_str(), "", _nBins, _min, _max, 100, 0., 1. );
  fullHist->SetXTitle( "#Lambda [TeV]" );
  fullHist->SetYTitle( _name.c_str() );
  fullHist->Draw();
  exclusion->Draw();
  fullHist->Draw( "ASAME" );

  TCanvas * zoomCanvas = new TCanvas( (_name+"_zoom").c_str(), "", 800, 800 );
  zoomCanvas->cd();
  TH2D * zoomHist = new TH2D( (_name+"zoomHist").c_str(), "", _nBins, _min, _max, 100, 0.04, 0.1 );
  zoomHist->SetXTitle( "#Lambda [TeV]" );
  zoomHist->SetYTitle( _name.c_str() );
  zoomHist->Draw();
  exclusion->Draw();
  zoomHist->Draw( "ASAME" );

}




ostream& operator<<( ostream& out, const CertaintyLevel& cl ) {

  out << cl._name << " observed: " << cl.observed() << endl;
  out << cl._name << " expected: " << cl.expected() << endl;

  return out;

}


double CertaintyLevel::observed() const {

  return findX( 0.05, _scales, _obs );

}

double CertaintyLevel::expected() const {
  return findX( 0.05, _scales, _exp );
}

double findX( const double& yVal, const vector<double>& x, const vector<double>& y ) {

  assert( x.size() == y.size() );
  double xHi = 0.;
  double xLo = 0.;
  double yHi = 0.;
  double yLo = 0.;
  
  for ( int i = 0; i < x.size(); ++i ) {
    if ( y.at(i) > yVal ) {
      xHi = x.at(i);
      xLo = x.at(i-1);
      yHi = y.at(i);
      yLo = y.at(i-1);
      break;
    }
  }

  // equation of a line ... solve for y
  double m = (yHi - yLo) / (xHi - xLo);

  return (yVal - yLo) / m + xLo;

}
