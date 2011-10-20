#include <iostream>
#include <vector>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TFitterMinuit.h>


#include "Likelihood.hpp"
#include "PDF.hpp"

using namespace std;


TH1 * CopyRange( TH1*, int, int );

int main() {

  cout << "hello frank!\n";

  TFile * pdfFile  = TFile::Open("~/docs/comp/analysis/vanilla.root");
  TFile * dataFile = TFile::Open("~/docs/comp/analysis/data.root");

  TH2 * pdfHist = (TH2*)pdfFile-> Get("PDF-2000-m_{jj}-7000 GeV");
  TH1 * fullHist= (TH1*)dataFile->Get("Chi_2000-to-7000all");
  TH1 * dataHist= CopyRange( fullHist, 1, 11 );
  
  TFitterMinuit * fitter = new TFitterMinuit();

  PDF pdf( pdfHist );
  vector<double> vec( 1, 1./(7000.*7000.) );
  cout << " pdf( 2, 5 | 0 ) =  " << pdf( 2., 5, vec ) << endl;

  Likelihood l( dataHist, &pdf );
  cout << " -2lnL( data | 0 ) = " << l( vec ) << endl;

  LikelihoodRatio launda( fitter, dataHist, &pdf );
  cout << " -2lnlaunda(0) = " << launda( vec ) << endl;

  return 0;

}


TH1 * CopyRange( TH1* h, int min, int max ) {
  vector<double> lowEdges;
  vector<double> content;
  vector<double> errors;

  for ( int i = min; i <= max; ++i ) {
    lowEdges.push_back( h->GetBinLowEdge( i ) );
    content. push_back( h->GetBinContent( i ) );
    errors.  push_back ( h->GetBinError( i ));
  }

  if ( max != h->GetNbinsX() ) lowEdges.push_back( h->GetBinLowEdge( max+1 ) );
  TH1* result = new TH1D( h->GetName(), "", lowEdges.size()-1, &lowEdges[0] );

  for ( int i = 0; i <= content.size(); ++i ) {
    result->SetBinContent( i+1, content[i] );
    result->SetBinError  ( i+1, errors[i]  ); 
  }

  return result;
  
}
