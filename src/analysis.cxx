#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

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

  TFile * pdfFile = TFile::Open( "~/docs/comp/analysis/kFactor.root" );
  TFile * dataFile = TFile::Open( "~/docs/comp/analysis/data.root" );

  TH1 * dataHist = (TH1*) dataFile->Get( "Chi_2000-to-7000all" );

  double nData = dataHist->Integral();
  cout << "nData = " << nData << endl;

  PDF pdf( pdfFile, dataHist->Integral() );
  pdf.useFit();
  PDFMonitor pdfMon;

  pdf.accept( pdfMon );

  PseudoExperimentFactory peFactory( &pdf, dataHist );
  vector< PseudoExperiment* > somePEs;
  vector<PseudoExperiment*> morePEs = peFactory.build( 0., 1.e4 );
  somePEs.insert( somePEs.end(), morePEs.begin(), morePEs.end() );
  vector< PseudoExperiment* > pValPEs;
  pValPEs.insert( pValPEs.end(), morePEs.begin(), morePEs.end() );

  morePEs = peFactory.build( 1 / pow( double( 2. ), 2 ), 1.e4 );
  somePEs.insert( somePEs.end(), morePEs.begin(), morePEs.end() );

  TestStatMonitor tm( "figures/Likelihood/", ".png" );
  foreach( PseudoExperiment* pe, somePEs )
  {
    Likelihood_FCN l( pe, &pdf, 1 / pow( double( 2. ), 2 ) );
    LikelihoodRatio_FCN launda( pe, &pdf, 1 / pow( double( 2. ), 2 ) );
    l.Minimize();

    l.accept( tm );
    launda.accept( tm );
  }

  tm.finalize();

  LikelihoodRatio_FCN lambda( dataHist, &pdf, 0. );
  PValueTest pv0( 0., lambda, pValPEs );

  TProfile * dataMinus2LogL = MapMinus2LogLikelihood( dataHist, pdf );
  TProfile * dataMinus2LogLambda = MapMinus2LogLikelihoodRatio( dataHist, pdf, 0. );
  TCanvas datac( "datac", "", 1000, 500 );
  datac.Divide( 2, 1 );
  datac.cd( 1 );
  dataMinus2LogL->Draw();
  datac.cd( 2 );
  dataMinus2LogLambda->Draw();
  datac.Print( "figures/dataMinus2LogL.png" );

  vector< PseudoExperiment* > someMorePEs = peFactory.build( 2., 1.e2 );
  TH1 * peHist = someMorePEs.at( 2 );
  TProfile * peMinus2LogL = MapMinus2LogLikelihood( peHist, pdf );
  TProfile * peMinus2LogLambda = MapMinus2LogLikelihoodRatio( peHist, pdf, 2. );
  TCanvas pec( "pec", "", 1000, 500 );
  pec.Divide( 2, 1 );
  pec.cd( 1 );
  peMinus2LogL->Draw();
  pec.cd( 2 );
  peMinus2LogLambda->Draw();
  pec.Print( "figures/peMinus2LogL.png" );

  peHist = somePEs.at( 2 );
  peMinus2LogL = MapMinus2LogLikelihood( peHist, pdf );
  peMinus2LogLambda = MapMinus2LogLikelihoodRatio( peHist, pdf, 0.25 );
  TCanvas pec2( "pec2", "", 1000, 500 );
  pec2.Divide( 2, 1 );
  pec2.cd( 1 );
  peMinus2LogL->Draw();
  pec2.cd( 2 );
  peMinus2LogLambda->Draw();
  pec2.Print( "figures/peMinus2LogL.png" );

  TCanvas * dataCanvas = new TCanvas( "dataCanvas", "", 500, 500 );
  dataCanvas->cd();
  dataHist->Draw();

  cout << " * pvalue = " << pv0( dataHist ) << endl;

  theApp.Run( kTRUE );
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
  TProfile* result = new TProfile( ( name + "_Minus2LogL" ).c_str(), "", 2000, 0, 4, 0., 1e6 );
  result->SetXTitle( "#alpha=1/#Lambda^{2}" );
  result->SetYTitle( ( "-2lnL(" + name + "|#alpha)" ).c_str() );

  vector< double > vec( 1, 0. );
  int nBins = result->GetNbinsX();
  double xMax = result->GetXaxis()->GetBinUpEdge( nBins );
  double delta = xMax / ( 2 * nBins );
  while ( vec.at( 0 ) < xMax ) {
    result->Fill( vec.at( 0 ), l( vec ) );
    vec.at( 0 ) += delta;
  }

  return result;

}

TProfile * MapMinus2LogLikelihoodRatio( TH1* exp, PDF& pdf, double alpha ) {

  LikelihoodRatio_FCN l( exp, &pdf, alpha );

  double max = 0.;
  double min = DBL_MAX;
  double maxAlpha = -1;
  double minAlpha = -1;

  string name = exp->GetName();
  TProfile* result = new TProfile( ( name + "_Minus2LogL" ).c_str(), "", 2000, 0, 4, 0., 1e6 );
  result->SetXTitle( "#alpha=1/#Lambda^{2}" );
  result->SetYTitle( ( name + " #lambda(#alpha)" ).c_str() );

  vector< double > vec( 1, 0. );
  int nBins = result->GetNbinsX();
  double xMax = result->GetXaxis()->GetBinUpEdge( nBins );
  double delta = xMax / ( 2 * nBins );
  while ( vec.at( 0 ) < xMax ) {
    result->Fill( vec.at( 0 ), l( vec ) );
    vec.at( 0 ) += delta;
  }

  return result;

}
