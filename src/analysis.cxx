#include <iostream>
#include <vector>
#include <cmath>
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

int main( int argc, char* argv[] ) {

  TApplication theApp( "Analysis", &argc, argv );
  SetAtlasStyle();
  // for GUI;
  TGClient windowClient;
  const TGWindow * rootWindow = windowClient.GetRoot();

  ControlFrame * control = new ControlFrame( rootWindow, 350, 80 );

  TFile * pdfFile = TFile::Open( "~/docs/comp/analysis/kFactor.root" );
  TFile * dataFile = TFile::Open( "~/docs/comp/analysis/data.root" );

  TGraph2D * pdfHist = (TGraph2D*) pdfFile->Get( "PDF-2000-m_{jj}-7000GeV" );
  //pdfHist->Smooth();
  TH1 * fullHist = (TH1*) dataFile->Get( "Chi_2000-to-7000all" );
  TH1 * dataHist = CopyRange( fullHist, 1, 11 );

  TCanvas * c = new TCanvas( "c", "", 500, 500 ); c->cd();
  dataHist->Draw();

  PDF pdf( pdfHist );
  PDFMonitor pdfMon;

  pdf.accept( pdfMon );

  PseudoExperimentFactory peFactory( &pdf, dataHist );
  vector< PseudoExperiment* > somePEs;
  //vector<PseudoExperiment*> morePEs = peFactory.build( 0., 1.e4 );
  //somePEs.insert( somePEs.end(), morePEs.begin(), morePEs.end() );

  vector< PseudoExperiment* > morePEs = peFactory.build( 1 / pow( double( 2. ), 2 ), 1.e3 );
  somePEs.insert( somePEs.end(), morePEs.begin(), morePEs.end() );

  TH1 * peHist = somePEs.at( 2 );

  vector< double > vec( 1, 0. );

  TestStatMonitor tm( "figures/Likelihood/", ".png" );
  foreach( PseudoExperiment* pe, somePEs ) {
    Likelihood_FCN l( pe, &pdf, 1 / pow( double( 2. ), 2 ) );
    LikelihoodRatio_FCN launda( pe, &pdf, 1 / pow( double( 2. ), 2 ) );

    l.accept( tm );
    launda.accept( tm );

  }

  tm.finalize();

  Likelihood_FCN l( dataHist, &pdf );
  LikelihoodRatio_FCN launda( dataHist, &pdf );

//  PValueTest pv0( 0., launda, somePEs );
//  cout << " * pvalue = " << pv0( dataHist ) << endl;

  double max = 0.;
  double min = 4e150;
  double maxAlpha = -1;
  double minAlpha = -1;

  TProfile dataMinus2LogL( "dataMinus2LogL", "", 5000, 0, 4, 0., 1e6 );
  dataMinus2LogL.SetXTitle( "#alpha=1/#Lambda^{2}" );
  dataMinus2LogL.SetYTitle( "-2lnL(data|#alpha)" );

  double delta = 4./(2*5000.);
  while ( vec.at( 0 ) < 4. ) {
    dataMinus2LogL.Fill( vec.at( 0 ), l( vec ) );
    if ( max - min < 1. ) {
      if ( max < l( vec ) ) {
        max = l( vec );
        maxAlpha = vec.at( 0 );
      }
      if ( min > l( vec ) ) {
        min = l( vec );
        minAlpha = vec.at( 0 );
      }
    }
    vec.at( 0 ) += delta;
  }

  cout << "max in min spread:\n" << " -2lnL(" << minAlpha << ") = " << min << endl << " -2lnL(" << maxAlpha << ") = "
       << max << endl << " |maxAlpha - minAlpha| = " << fabs( maxAlpha - minAlpha ) << endl << " |max - min| = "
       << fabs( max - min ) << endl;

  TCanvas datac( "datac", "", 500, 500 );
  datac.cd();
  //datac.SetLogy();
  dataMinus2LogL.Draw();
  datac.Print( "figures/dataMinus2LogL.png" );
  vec.at( 0 ) = 0;

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
