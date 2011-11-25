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

#include <TApplication.h>
#include <TCanvas.h>
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

  // Set up PDF for data
  PDF pdf( pdfFile, dataHist->Integral() );
  pdf.useFit();

  // Monitor data PDF
  PDFMonitor pdfMon;
  pdf.accept( pdfMon );

  int nPE = 1.e4;
  vector< double > alphas = list_of( 0. )( pow( 1. / 8., 4 ) )( pow( 1. / 7., 4 ) )( pow( 1. / 6., 4 ) )(
      pow( 1. / 5., 4 ) )( pow( 1. / 4., 4 ) )( pow( 1. / 3., 4 ) )( pow( 1. / 1.5, 4 ) )( pow( 1. / 1., 4 ) )(
      pow( 1. / .75, 4 ) )( pow( 1. / .5, 4 ) );
  PseudoExperimentFactory peFactory( &pdf, dataHist );
  typedef map< double, vector< PseudoExperiment* > > peMap_t;
  peMap_t pes;
  foreach( const double& alpha, alphas )
    pes[alpha] = peFactory.build( alpha, nPE );

  foreach( const peMap_t::value_type& x, pes )
  {
    const double alpha = x.first;
    TestStatMonitor tm( alpha, "figures/Likelihood/", ".png" );
    vector< LikelihoodRatio* > likelihoodRatios;
    likelihoodRatios.reserve( x.second.size() );
    foreach( PseudoExperiment* pe, x.second )
    {
      PDF * pePDF = new PDF( pdfFile, pe->Integral() );
      likelihoodRatios.push_back( new LikelihoodRatio( pe, pePDF, alpha ) );
    }
    foreach( LikelihoodRatio* lambda, likelihoodRatios )
    {
      lambda->accept( tm );
      lambda->denominator().accept( tm );
    }
    tm.finalize();
  }

  double alpha = alphas.front();
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

  vector< LikelihoodRatio* > likelihoodRatios;
  likelihoodRatios.reserve( nPE );
  foreach( PseudoExperiment* pe, pes[alpha] )
  {
    PDF * pePDF = new PDF( pdfFile, pe->Integral() );
    likelihoodRatios.push_back( new LikelihoodRatio( pe, pePDF, alpha ) );
  }
  PValueTest pv0( alpha, likelihoodRatios );
  LikelihoodRatio dataLikelihoodRatio( dataHist, new PDF( pdfFile, dataHist->Integral() ), alpha );
  double pValue = pv0( dataLikelihoodRatio );
  cout << " * pvalue( Lambda = " << 1. / pow( alpha, 0.25 ) << ") = " << pValue << endl;
  pv0.finalize();

  theApp.Run( kTRUE );

  pdfFile->Close();
  dataFile->Close();
  //for_each( likelihoodRatios.begin(), likelihoodRatios.end(), boost::checked_deleter< LikelihoodRatio >() );
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
