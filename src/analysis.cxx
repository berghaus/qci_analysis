#include <iostream>
#include <vector>
#include <cmath>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


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

#include "PDFMonitor.hpp"
#include "TestStatMonitor.hpp"

using namespace std;


TH1 * CopyRange( TH1*, int, int );

int main() {

  //  cout << "hello frank!\n";
  SetAtlasStyle();

  TFile * pdfFile  = TFile::Open("~/docs/comp/analysis/vanilla.root");
  TFile * dataFile = TFile::Open("~/docs/comp/analysis/data.root");

  TH2 * pdfHist = (TH2*)pdfFile-> Get("PDF-2000-m_{jj}-7000 GeV");
  pdfHist->Smooth();
  TH1 * fullHist= (TH1*)dataFile->Get("Chi_2000-to-7000all");
  TH1 * dataHist= CopyRange( fullHist, 1, 11 );
  
  PDF pdf( pdfHist );
  PDFMonitor pdfMon;

  pdf.accept( pdfMon );

  PseudoExperimentFactory peFactory( &pdf, dataHist );
  vector<PseudoExperiment*> somePEs = peFactory.build( 0., 1.e4 );
  
  TH1 * peHist = somePEs.at(2);

  vector<double> vec( 1, 0. );

  TestStatMonitor tm("figures/Likelihood/",".png" );
  foreach( PseudoExperiment* pe, somePEs ) {
    Likelihood_FCN      l( pe, &pdf );
    LikelihoodRatio_FCN launda( pe, &pdf );

    l.     accept( tm );
    launda.accept( tm );
    
  }
  Likelihood_FCN  l( dataHist, &pdf );
  LikelihoodRatio_FCN launda( dataHist, &pdf );

  // PValueTest pv0( 0., launda, somePEs );
  // cout << " * pvalue = " <<  pv0( dataHist ) << endl;


  // double max = 0.;
  // double min = 4e150;
  // double maxAlpha = -1 ;
  // double minAlpha = -1 ;

  // TProfile dataMinus2LogL("dataMinus2LogL","",500,0,1e-7,50,1e2);
  // dataMinus2LogL.SetXTitle("#alpha=1/#Lambda^{2}");
  // dataMinus2LogL.SetYTitle("-2lnL(data|#alpha)");

  // double delta = 1.e-11;
  // while ( vec.at(0) < 1.e-7 ) {
  //   if ( isinf( l(vec) ) || isnan( l(vec) )  ) {
  //     cout << "PDF broke." << endl; 
  //     break;
  //   }
  //   dataMinus2LogL.Fill( vec.at(0), l( vec ) );
  //   if ( max - min < 1. ) {
  //     if ( max < l( vec ) ) { max = l( vec ); maxAlpha = vec.at(0); }
  //     if ( min > l( vec ) ) { min = l( vec ); minAlpha = vec.at(0); }
  //   }
  //   vec.at(0) += delta;
  // }

  // cout << "max in min spread:\n" 
  //      << " -2lnL(" << minAlpha << ") = " << min << endl
  //      << " -2lnL(" << maxAlpha << ") = " << max << endl
  //      << " |maxAlpha - minAlpha| = " << fabs(maxAlpha-minAlpha) << endl
  //      << " |max - min| = " << fabs(max-min) << endl;

  // TCanvas c("c","",800,600); c.cd(); c.SetLogy();
  // dataMinus2LogL.Draw();
  // c.Print("figures/dataMinus2LogL.png");
  // vec.at(0) = 0;

  return 0;

}


TH1 * CopyRange( TH1* h, int min, int max ) {
  vector<double> lowEdges;
  vector<double> content;
  vector<double> errors;

  for ( int i = min; i <= max; ++i ) {
    lowEdges.push_back( h->GetBinLowEdge( i ) );
    content. push_back( h->GetBinContent( i ) );
    errors.  push_back( h->GetBinError( i ) );
  }

  if ( max != h->GetNbinsX() ) lowEdges.push_back( h->GetBinLowEdge( max+1 ) );
  TH1* result = new TH1D( h->GetName(), "", lowEdges.size()-1, &lowEdges[0] );

  for ( int i = 0; i <= content.size(); ++i ) {
    result->SetBinContent( i+1, content[i] );
    result->SetBinError  ( i+1, errors[i]  ); 
  }

  return result;
  
}
