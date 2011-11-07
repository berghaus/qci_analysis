#include "PDFMonitor.hpp"
#include <vector>
#include <boost/format.hpp>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

#include "PDF.hpp"
using namespace std;
using boost::format;

PDFMonitor::PDFMonitor() :
    _folder( "figures/PDF/" ),
    _ext( ".png" ) {
}

PDFMonitor::PDFMonitor( const string& folder, const string ext ) :
    _folder( folder ),
    _ext( ext ) {
}

PDFMonitor::~PDFMonitor() {
}

void PDFMonitor::monitor( PDF& pdf ) {

  // make some histograms
  TCanvas *pdfCanvas = new TCanvas( "PDFMonCanvas", "", 800, 600 );
  pdfCanvas->cd();
  pdfCanvas->Divide( 4, 3 );

  for( int xBin = 1; xBin <= pdf.hist()->GetNbinsX(); ++xBin ) {

    double chi = pdf.hist()->GetXaxis()->GetBinCenter( xBin );
    if ( chi < 0. ) continue;
    if ( 30. < chi ) break;

    string projName = string( pdf.hist()->GetName() ) + str( format( "-chi-%2.1f" ) % chi );
    TH1* proj = pdf.hist()->ProjectionY( projName.c_str(), xBin, xBin );

    // interpolation vector
    vector< double > alphas( 1, 0. );
    vector< double > interpol( 1, pdf.interpolate( chi, alphas.back() ) );

    double delta = 1.e-3;
    while ( alphas.back() < 4. ) {
      alphas.push_back( alphas.back() + delta );
      interpol.push_back( pdf.interpolate( chi, alphas.back() ) );
    }
    TGraph * interpolation = new TGraph( alphas.size(), &alphas[0], &interpol[0] );
    interpolation->SetLineColor( kBlue );

    pdfCanvas->cd(xBin);
    proj->Draw( "EX0" );
    interpolation->Draw( "SAMEL" );
  }
  pdfCanvas->Print( ( _folder + "PDFs" + _ext ).c_str() );

}

