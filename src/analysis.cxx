#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <limits>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/checked_delete.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/math/distributions/normal.hpp>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH2.h>
#include <TFile.h>

#include "Likelihood.hpp"
#include "PDF.hpp"
#include "PseudoExperiment.hpp"
#include "PValueTest.hpp"
#include "AtlasStyle.hpp"
#include "ControlFrame.hpp"

#include "PDFMonitor.hpp"
#include "TestStatMonitor.hpp"

using namespace std;
using namespace boost::assign;
using namespace boost::math;

template< class T >
T quantile( const vector<T>&, const double& );
TH1 * CopyRange( TH1*, int, int );
TProfile * MapMinus2LogLikelihood( TH1*, PDF*, double );
TProfile * MapMinus2LogLikelihoodRatio( TH1*, PDF*, double );

int main( int argc, char* argv[] ) {

  // for GUI;
  TApplication theApp( "Analysis", &argc, argv );
  SetAtlasStyle();
  TGClient windowClient;
  const TGWindow * rootWindow = windowClient.GetRoot();
  ControlFrame * control = new ControlFrame( rootWindow, 350, 80 );

  TFile * pdfFile = TFile::Open( "~/docs/comp/analysis/kFactor.root", "READ" );
  TFile * dataFile = TFile::Open( "~/docs/comp/analysis/data.root", "READ" );
  TH1 * dataHist = (TH1*) dataFile->Get( "Chi_2000-to-7000all" );

  normal myGaus;

  // Set up PDF for data
  PDF pdf( pdfFile, dataHist->Integral() );
  pdf.useFit();

  LikelihoodRatio dataLikelihoodRatio( dataHist, new PDF( pdfFile, dataHist->Integral() ), 0. );

  // Monitor data PDF
  PDFMonitor pdfMon;
  pdf.accept( pdfMon );

  PseudoExperimentFactory peFactory( &pdf, dataHist );

  // TODO: give these as command line controls
  int nPE = 10000;
  double nBinsScale = 50;
  double minScale   = 4.;
  double maxScale   = 8.;
  double deltaScale = ( maxScale - minScale ) / nBinsScale;
  typedef map< double, vector< PseudoExperiment* > > peMap_t;
  peMap_t peMap;
  vector<PValueTest*> pValueTests;
  vector<double> scales;
  vector<double> observed;
  vector<double> expected;
  vector<double> e02;
  vector<double> e98;
  vector<double> e16;
  vector<double> e84;
  vector<double> ex;

  vector< PseudoExperiment* > errorBandPEs = peFactory.build( 0., nPE );
  vector< LikelihoodRatio* > errorBandLRs;
  errorBandLRs.reserve( errorBandLRs.size() );
  foreach( PseudoExperiment* pe, errorBandPEs ) {
    PDF * pePDF = new PDF( pdfFile, pe->Integral() );
    errorBandLRs.push_back( new LikelihoodRatio( pe, pePDF, 0. ) );
  }

  PValueTest pValueTest( 0., errorBandLRs );
  double pValue = pValueTest( dataLikelihoodRatio );
  cout << " * pvalue( Lambda = " << 0. << ") = " << pValue << endl;  

  pValueTest.finalize();

  theApp.Run( kTRUE );
  return 0;

  for ( double scale = minScale; scale < maxScale; scale += deltaScale ) {
    double alpha = pow( scale, -4 );
    vector< PseudoExperiment* > pes = peFactory.build( alpha, nPE );
    peMap[alpha] = pes;

    vector< LikelihoodRatio* > likelihoodRatios;
    likelihoodRatios.reserve( pes.size() );
    foreach( PseudoExperiment* pe, pes )
    {
      PDF * pePDF = new PDF( pdfFile, pe->Integral() );
      likelihoodRatios.push_back( new LikelihoodRatio( pe, pePDF, alpha ) );
    }
    PValueTest * pValueTest = new PValueTest( alpha, likelihoodRatios );
    PValueTest& pv = *pValueTest;
    pValueTests.push_back( pValueTest );
    double pValue = pv( dataLikelihoodRatio );
    cout << " * pvalue( Lambda = " << scale << ") = " << pValue << endl;

    vector< double > ebPValues; ebPValues.reserve( nPE );
    foreach( LikelihoodRatio* l, errorBandLRs ) ebPValues.push_back( pv( *l ) );
    sort( ebPValues.begin(), ebPValues.end() );

    double sigmaN2 = quantile( ebPValues, cdf( myGaus, -2. ) );
    double sigmaN1 = quantile( ebPValues, cdf( myGaus, -1. ) );
    double exp = quantile( ebPValues, cdf( myGaus, 0. ) );
    double sigmaP1 = quantile( ebPValues, cdf( myGaus, 1. ) );
    double sigmaP2 = quantile( ebPValues, cdf( myGaus, 2. ) );

    scales.push_back( scale );
    observed.push_back( pValue );
    expected.push_back( exp );
    e02.push_back( fabs( exp - sigmaN2 ) );
    e98.push_back( fabs( exp - sigmaP2 ) );
    e16.push_back( fabs( exp - sigmaN1 ) );
    e84.push_back( fabs( exp - sigmaP1 ) );
    ex.push_back( 0. );
    for_each( likelihoodRatios.begin(), likelihoodRatios.end(), boost::checked_deleter< LikelihoodRatio >() );
  }

  TGraph observedGraph( scales.size(), &scales[0], &observed[0] );
  observedGraph.SetLineWidth( 2 );
  observedGraph.SetLineColor( kRed );

  TGraphAsymmErrors expected1Sigma( nBinsScale, &scales[0], &expected[0], &ex[0], &ex[0], &e16[0], &e84[0] );
  expected1Sigma.SetFillStyle( 1001 );
  expected1Sigma.SetFillColor( kBlue );
  expected1Sigma.SetLineStyle( 2 );
  expected1Sigma.SetLineWidth( 2 );

  TGraphAsymmErrors expected2Sigma( nBinsScale, &scales[0], &expected[0], &ex[0], &ex[0], &e02[0], &e98[0] );
  expected2Sigma.SetFillStyle( 1001 );
  expected2Sigma.SetFillColor( kYellow );
  expected2Sigma.SetLineColor( kYellow );
  expected2Sigma.SetLineWidth( 2 );

  TMultiGraph exclusion;
  exclusion.Add( &expected2Sigma, "3" );
  exclusion.Add( &expected1Sigma, "3" );
  exclusion.Add( &expected1Sigma, "LX" );
  exclusion.Add( &observedGraph, "L" );

  //  PseudoExperiment * pe = pes[alpha].front();
  //  PDF * pePDF = new PDF( pdfFile, pe->Integral() );
  //  TProfile * peMinus2LogL = MapMinus2LogLikelihood( pe, pePDF, alpha );
  //  TProfile * peMinus2LogLambda = MapMinus2LogLikelihoodRatio( pe, pePDF, alpha );
  //  TCanvas pec( "pec", "", 1000, 500 );
  //  pec.Divide( 2, 1 );
  //  pec.cd( 1 );
  //  peMinus2LogL->Draw();
  //  pec.cd( 2 );
  //  peMinus2LogLambda->Draw();
  //  pec.Print( "figures/peMinus2LogL.png" );

  TCanvas obsCanvas( "observedCanvas", "", 800, 800 );
  obsCanvas.cd();
  exclusion.Draw("a");
  

  theApp.Run( kTRUE );

  pdfFile->Close();
  dataFile->Close();
  return 0;

}

TH1 * CopyRange( TH1* h, int min, int max ) {
  vector< double > lowEdges;
  vector< double > content;
  vector< double > errors;

  for( int i = min; i <= max; ++i ) {
    lowEdges.push_back( h->GetBinLowEdge( i ) );
    content.push_back( h->GetBinContent( i ) );
    errors.push_back( h->GetBinError( i ) );
  }

  if ( max != h->GetNbinsX() ) lowEdges.push_back( h->GetBinLowEdge( max + 1 ) );
  TH1* result = new TH1D( h->GetName(), "", lowEdges.size() - 1, &lowEdges[0] );

  for( int i = 0; i <= content.size(); ++i ) {
    result->SetBinContent( i + 1, content[i] );
    result->SetBinError( i + 1, errors[i] );
  }

  return result;

}

TProfile * MapMinus2LogLikelihood( TH1* exp, PDF* pdf, double alpha ) {

  Likelihood_FCN l( exp, pdf );

  double offset = 0.1;
  double min = alpha / 16 + offset;
  double max = alpha * 16 + offset;
  double nBins = 1000;
  vector< double > alphaBins;
  alphaBins.reserve( nBins );
  for( int iBin = 0; iBin <= nBins; ++iBin ) {
    double binEdge = pow( max, double( iBin ) / double( nBins ) )
                     * pow( min, double( nBins - iBin ) / double( nBins ) );
    alphaBins.push_back( binEdge );
  }

  string name = exp->GetName();
  TProfile* result = new TProfile( ( name + "_Minus2LogL" ).c_str(), "", alphaBins.size() - 1, &alphaBins[0], -10.,
                                   250. );
  result->SetXTitle( "#alpha + 0.1 [TeV^{-4}]" );
  result->SetYTitle( ( " L(  " + name + " | #alpha )" ).c_str() );

  for( int i = 0; i < alphaBins.size() - 1; ++i ) {
    double a = alphaBins[i];
    double xMax = alphaBins[i + 1];
    double delta = ( xMax - a ) / 5.;
    while ( a < xMax ) {
      vector< double > vec( 1, a - offset );
      result->Fill( a, l( vec ) );
      a += delta;
    }
  }

  return result;

}

TProfile * MapMinus2LogLikelihoodRatio( TH1* exp, PDF* pdf, double alpha ) {

  LikelihoodRatio l( exp, pdf, alpha );

  double offset = 0.1;
  double min = alpha / 16 + offset;
  double max = alpha * 16 + offset;
  double nBins = 1000;
  vector< double > alphaBins;
  alphaBins.reserve( nBins );
  for( int iBin = 0; iBin <= nBins; ++iBin ) {
    double binEdge = pow( max, double( iBin ) / double( nBins ) )
                     * pow( min, double( nBins - iBin ) / double( nBins ) );
    alphaBins.push_back( binEdge );
  }

  string name = exp->GetName();
  TProfile* result = new TProfile( ( name + "_Minus2LogLRatio" ).c_str(), "", alphaBins.size() - 1, &alphaBins[0], -10.,
                                   250. );
  result->SetXTitle( "#alpha + 0.1 [TeV^{-4}]" );
  result->SetYTitle( ( name + " #lambda( #alpha )" ).c_str() );

  for( int i = 0; i < alphaBins.size() - 1; ++i ) {
    double a = alphaBins[i];
    double xMax = alphaBins[i + 1];
    double delta = ( xMax - a ) / 5.;
    while ( a < xMax ) {
      vector< double > vec( 1, a - offset );
      result->Fill( a, l( vec ) );
      a += delta;
    }
  }

  return result;

}


template< class T >
T quantile( const vector<T>& v, const double& q ) {
  typename vector<T>::size_type index = v.size() * q;
  return v.at( index );
}

