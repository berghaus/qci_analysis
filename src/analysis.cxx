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

TH1 * CopyRange( TH1*, int, int );
TProfile * MapMinus2LogLikelihood( TH1*, PDF& );
TProfile * MapMinus2LogLikelihoodRatio( TH1*, PDF&, double );

int main( int argc, char* argv[] ) {

  TApplication theApp( "Analysis", &argc, argv );
  SetAtlasStyle();
  // for GUI;
  TGClient windowClient;
  const TGWindow * rootWindow = windowClient.GetRoot();

  ControlFrame * control = new ControlFrame( rootWindow, 350, 80 );

  TFile * pdfFile = TFile::Open( "~/docs/comp/analysis/kFactor.root", "READ" );
  TFile * dataFile = TFile::Open( "~/docs/comp/analysis/data.root", "READ" );

  TH1 * dataHist = (TH1*) dataFile->Get( "Chi_2000-to-7000all" );

  PDF pdf( pdfFile, dataHist->Integral() );
  pdf.useFit();
  PDFMonitor pdfMon;

  pdf.accept( pdfMon );

  double alpha = 0; // TeV
  PseudoExperimentFactory peFactory( &pdf, dataHist );
  vector< PseudoExperiment* > somePEs;
  vector< PseudoExperiment* > morePEs = peFactory.build( alpha, 1.e2 );
  somePEs.insert( somePEs.end(), morePEs.begin(), morePEs.end() );
  vector< PseudoExperiment* > pValPEs = morePEs;

  TestStatMonitor tm( "figures/Likelihood/", ".png" );

  vector< LikelihoodRatio* > likelihoodRatios;
  likelihoodRatios.reserve( somePEs.size() );
  foreach( PseudoExperiment* pe, somePEs )
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

//  TProfile * dataMinus2LogL = MapMinus2LogLikelihood( dataHist, pdf );
//  TProfile * dataMinus2LogLambda = MapMinus2LogLikelihoodRatio( dataHist, pdf, 0. );
//  TCanvas datac( "datac", "", 1000, 500 );
//  datac.Divide( 2, 1 );
//  datac.cd( 1 );
//  dataMinus2LogL->Draw();
//  datac.cd( 2 );
//  dataMinus2LogLambda->Draw();
//  datac.Print( "figures/dataMinus2LogL.png" );

  PValueTest pv0( alpha, likelihoodRatios );
  LikelihoodRatio dataLikelihoodRatio( dataHist, new PDF( pdfFile, dataHist->Integral() ), alpha );
  double pValue = pv0( dataLikelihoodRatio );
  cout << " * pvalue( Lambda = " << 1. / pow( alpha, 0.25 ) << ") = " << pValue << endl;
  pv0.finalize();

  theApp.Run( kTRUE );

  pdfFile->Close();
  dataFile->Close();
  for_each( likelihoodRatios.begin(), likelihoodRatios.end(), boost::checked_deleter< LikelihoodRatio >() );
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

TProfile * MapMinus2LogLikelihood( TH1* exp, PDF& pdf ) {

  Likelihood_FCN l( exp, &pdf );

  double max = 0.;
  double min = DBL_MAX;
  double maxAlpha = -1;
  double minAlpha = -1;

  string name = exp->GetName();
  TProfile* result = new TProfile( ( name + "_Minus2LogL" ).c_str(), "", 2000, 0.4, 20, 0., 1e6 );
  result->SetXTitle( "#Lambda [TeV]" );
  result->SetYTitle( ( "-2lnL( " + name + " | #Lambda )" ).c_str() );

  double scale = result->GetXaxis()->GetBinLowEdge( 1 );
  int nBins = result->GetNbinsX();
  double xMax = result->GetXaxis()->GetBinUpEdge( nBins );
  double delta = xMax / ( 2 * nBins );
  while ( scale < xMax ) {
    vector< double > vec( 1, 1. / pow( scale, 4 ) );
    result->Fill( scale, l( vec ) );
    scale += delta;
  }

  return result;

}

TProfile * MapMinus2LogLikelihoodRatio( TH1* exp, PDF& pdf, double alpha ) {

  LikelihoodRatio l( exp, &pdf, alpha );

  double max = 0.;
  double min = DBL_MAX;
  double maxAlpha = -1;
  double minAlpha = -1;

  string name = exp->GetName();
  TProfile* result = new TProfile( ( name + "_Minus2LogL" ).c_str(), "", 2000, 0.4, 20, 0., 1e6 );
  result->SetXTitle( "#Lambda [TeV]" );
  result->SetYTitle( ( name + " #lambda( #Lambda )" ).c_str() );

  double scale = result->GetXaxis()->GetBinLowEdge( 1 );
  int nBins = result->GetNbinsX();
  double xMax = result->GetXaxis()->GetBinUpEdge( nBins );
  double delta = xMax / ( 2 * nBins );
  while ( scale < xMax ) {
    vector< double > vec( 1, 1. / pow( scale, 4 ) );
    result->Fill( scale, l( vec ) );
    scale += delta;
  }

  return result;

}
