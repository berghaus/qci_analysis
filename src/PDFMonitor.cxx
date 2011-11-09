#include "PDFMonitor.hpp"
#include <map>
#include <string>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

#include "PDF.hpp"
using namespace std;
#define foreach BOOST_FOREACH
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
  int nPad = 1;

  typedef map< double, TGraphErrors* > chiGraphMap_t;
  chiGraphMap_t graphs = pdf.eventCounts();
  foreach( chiGraphMap_t::value_type ec, graphs )
  {
    double chi = ec.first;
    TGraphErrors *graph = ec.second;
    string title = str( format(  "#chi = %2.1f" ) % chi );

    pdfCanvas->cd( nPad )->SetLogx();
    graph->SetTitle( title.c_str() );
    graph->Draw("AP");

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
    interpolation->Draw( "SAMEL" );
    ++nPad;
  }

  pdfCanvas->Print( ( _folder + "PDFs" + _ext ).c_str() );

}

