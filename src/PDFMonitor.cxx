#include "PDFMonitor.hpp"
#include <algorithm>
#include <functional>
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
    _ext( ".png" ),
    _pdfCanvas( "PDFMonCanvas", "", 800, 600 ),
    _interpolCanvas( "InterpolCanvas", "", 800, 600 ),
    _fitResultCanvas( "FitResultCanvas", "", 800, 600 ) {
  _pdfCanvas.Divide( 4, 3 );
  _interpolCanvas.Divide( 4, 3 );
  _fitResultCanvas.Divide( 4, 3 );
}

PDFMonitor::PDFMonitor( const string& folder, const string ext ) :
    _folder( folder ),
    _ext( ext ),
    _pdfCanvas( "PDFMonCanvas", "", 800, 600 ),
    _interpolCanvas( "InterpolCanvas", "", 800, 600 ),
    _fitResultCanvas( "FitResultCanvas", "", 800, 600 ) {
  _pdfCanvas.Divide( 4, 3 );
  _interpolCanvas.Divide( 4, 3 );
  _fitResultCanvas.Divide( 4, 3 );
}

PDFMonitor::~PDFMonitor() {
  // delete _interpolations
  for_each( _interpolations.begin(), _interpolations.end(), bind2nd( mem_fun( &TGraph::Delete ), "" ) );
}

void PDFMonitor::monitor( PDF& pdf ) {

  // make some histograms
  int nPad = 1;

  typedef map< double, TGraphErrors* > chiGraphMap_t;
  chiGraphMap_t graphs = pdf.eventCounts();
  foreach( chiGraphMap_t::value_type ec, graphs )
  {
    double chi = ec.first;
    TGraphErrors *graph = ec.second;
    string title = str( format( "#chi = %2.1f" ) % chi );
    TVirtualPad * pad = _pdfCanvas.cd( nPad );
    TF1 * func = graph->GetFunction( "PDFFit" );
    func->SetLineColor( kRed );
    graph->SetTitle( title.c_str() );
    graph->Draw( "AP" );
    graph->GetXaxis()->SetTitle( "#alpha = 1/#Lambda^{4} [TeV^{-4}]" );
    graph->GetYaxis()->SetTitle( "n(#alpha)" );
    ++nPad;
  }

  nPad = 1;
  foreach( chiGraphMap_t::value_type ec, graphs )
  {
    double chi = ec.first;
    string title = str( format( "#chi = %2.1f" ) % chi );
    string fname = str( format( "FunChi%2.1f" ) % chi );
    TVirtualPad * pad = _interpolCanvas.cd( nPad );

    vector< double > scales( 1, 0.1 );
    vector< double > interpol( 1, pdf.interpolate( chi, 1. / pow( scales.back(), 4 ) ) );
    double delta = 0.1;
    while ( scales.back() < 20. ) {
      scales.push_back( scales.back() + delta );
      interpol.push_back( pdf.interpolate( chi, 1. / pow( scales.back(), 4 ) ) );
    }
    TGraph * interpolation = new TGraph( scales.size(), &scales[0], &interpol[0] );
    _interpolations.push_back( interpolation );
    interpolation->SetLineColor( kBlue );
    interpolation->SetTitle( title.c_str() );
    interpolation->GetXaxis()->SetTitle( "#Lambda [TeV]" );
    interpolation->GetYaxis()->SetTitle( "n(#Lambda)" );
    interpolation->Draw( "AL" );
    ++nPad;
  }

  nPad = 1;
  foreach( chiGraphMap_t::value_type ec, graphs )
  {
    double chi = ec.first;
    TGraphErrors *graph = ec.second;
    string title = str( format( "#chi = %2.1f" ) % chi );
    TVirtualPad * pad = _fitResultCanvas.cd( nPad );
    TF1 * func = graph->GetFunction( "PDFFit" );
    int n = graph->GetN();
    double * xArr = graph->GetX();
    double * yArr = graph->GetY();
    double * eyArr = graph->GetEY();
    vector< double > index;
    vector< double > offset;
    index.reserve( n );
    offset.reserve( n );
    for( int i = 0; i < n; ++i ) {
      index.push_back( i );
      double graphY = yArr[i];
      double funcY = func->Eval( xArr[i] );
      offset.push_back( ( graphY - funcY ) / eyArr[i] );
    }
    TGraph * fitResult = new TGraph( n, &index[0], &offset[0] );
    _fitResults.push_back( fitResult );
    fitResult->SetMarkerColor( kRed );
    fitResult->SetTitle( title.c_str() );
    fitResult->GetXaxis()->SetTitle( "Index of #alpha = #Lambda^{-4}" );
    fitResult->GetYaxis()->SetTitle( "( n^{MC}(#alpha) - n^{fit}(#alpha) ) / #sigma^{MC}" );
    fitResult->Draw( "AP" );
    ++nPad;
  }

  std::vector< TGraph* > _fitResults;

  _pdfCanvas.Print( ( _folder + "PDFs" + _ext ).c_str() );
  _interpolCanvas.Print( ( _folder + "Interpolation" + _ext ).c_str() );
  _fitResultCanvas.Print( ( _folder + "FitResults" + _ext ).c_str() );

}

