#include "PDFMonitor.hpp"
#include <cmath>
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
#include "TF1.h"

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
  TCanvas *interpolCanvas = new TCanvas( "InterpolCanvas", "", 800, 600 );
  pdfCanvas->Divide( 4, 3 );
  interpolCanvas->Divide( 4, 3 );

  pdfCanvas->cd();
  int nPad = 1;

  typedef map< double, TGraphErrors* > chiGraphMap_t;
  chiGraphMap_t graphs = pdf.eventCounts();
  foreach( chiGraphMap_t::value_type ec, graphs )
  {
    double chi = ec.first;
    TGraphErrors *graph = ec.second;
    string title = str( format( "#chi = %2.1f" ) % chi );
    TVirtualPad * pad = pdfCanvas->cd( nPad );
    TF1 * func = graph->GetFunction( "PDFFit" );
    func->SetLineColor( kRed );
    graph->SetTitle( title.c_str() );
    graph->Draw( "AP" );
    graph->GetXaxis()->SetTitle( "#alpha = 1/#Lambda^{4} [TeV^{-4}]" );
    graph->GetYaxis()->SetTitle( "n(#alpha)" );
    ++nPad;
  }

  interpolCanvas->cd();
  nPad = 1;
  foreach( chiGraphMap_t::value_type ec, graphs )
  {
    double chi = ec.first;
    string title = str( format( "#chi = %2.1f" ) % chi );
    string fname = str( format( "FunChi%2.1f" ) % chi );
    TVirtualPad * pad = interpolCanvas->cd( nPad );

    vector< double > scales( 1, 0.1 );
    vector< double > interpol( 1, pdf.interpolate( chi, 1./pow(scales.back(),4) ) );
    double delta = 0.1;
    while ( scales.back() < 20. ) {
      scales.push_back( scales.back() + delta );
      interpol.push_back( pdf.interpolate( chi, 1./pow(scales.back(),4) ) );
    }
    TGraph * interpolation = new TGraph( scales.size(), &scales[0], &interpol[0] );
    interpolation->SetLineColor( kBlue );
    interpolation->SetTitle( title.c_str() );
    interpolation->GetXaxis()->SetTitle( "#Lambda [TeV]" );
    interpolation->GetYaxis()->SetTitle( "n(#Lambda)" );
    interpolation->Draw( "AL" );
    ++nPad;
  }

  pdfCanvas->Print( ( _folder + "PDFs" + _ext ).c_str() );
  interpolCanvas->Print( ( _folder + "Interpolation" + _ext ).c_str() );

}

