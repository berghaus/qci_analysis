#include "PredictionMonitor.hpp"
#include <algorithm>
#include <functional>
#include <cmath>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/checked_delete.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/list_of.hpp>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"


#include "Prediction.hpp"
using namespace std;
#define foreach BOOST_FOREACH
using namespace boost;

template<class C>
string CastToString( C& n ){
  ostringstream ss;
  ss << fixed << setprecision(2) << n;
  return ss.str();
}


//_____________________________________________________________________________________________________________________
PredictionMonitor::PredictionMonitor() :
    _folder( "figures/PDF/" ),
    _ext( ".pdf" ),
    _pdfCanvas( "PDFMonCanvas", "", 1800, 2400 ),
    _interpolCanvas( "InterpolCanvas", "", 1800, 2400 ),
    _fitResultCanvas( "FitResultCanvas", "", 1800, 2400 ),
    _parameterCanvas( "ParameterCanvas", "", 1800, 2400 ) {
  _pdfCanvas.Divide( 3, 4, 0, 0 );
  _interpolCanvas.Divide( 3, 4, 0, 0 );
  _fitResultCanvas.Divide( 3, 4, 0, 0 );
}

//_____________________________________________________________________________________________________________________
PredictionMonitor::PredictionMonitor( const string& folder, const string ext ) :
    _folder( folder ),
    _ext( ext ),
    _pdfCanvas( "PDFMonCanvas", "", 1800, 2400 ),
    _interpolCanvas( "InterpolCanvas", "", 1800, 2400 ),
    _fitResultCanvas( "FitResultCanvas", "", 1800, 2400 ),
    _parameterCanvas( "ParameterCanvas", "", 1800, 2400 ) {
  _pdfCanvas.Divide( 3, 4, 0, 0 );
  _interpolCanvas.Divide( 3, 4, 0, 0 );
  _fitResultCanvas.SetBottomMargin(0.2);
  _fitResultCanvas.Divide( 3, 4, 0, 0 );
}

//_____________________________________________________________________________________________________________________
PredictionMonitor::~PredictionMonitor() {
  // delete _interpolations
  for_each( _interpolations.begin(), _interpolations.end(), checked_deleter< TGraph >() );
  for_each( _fitResults.begin(), _fitResults.end(), checked_deleter< TGraph >() );
}

//_____________________________________________________________________________________________________________________
void PredictionMonitor::monitor( Prediction& pdf ) {

  double mjj = 2000.; // need to re-implement this visitor
  // make some histograms
  int nPad = 1;
  typedef map< double, TGraphErrors* > chiGraphMap_t;
  chiGraphMap_t graphs = pdf.eventCounts(mjj);
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
    pdf.setMjj( mjj );
    while ( scales.back() < 20. ) {
      scales.push_back( scales.back() + delta );
      interpol.push_back( pdf.interpolate( chi, 1. / pow( scales.back(), 4 ) ) );
    }
    TGraph * interpolation = new TGraph( scales.size(), &scales[0], &interpol[0] );
    _interpolations.push_back( interpolation );
    interpolation->SetLineColor( kBlue );
    interpolation->SetTitle( title.c_str() );
    interpolation->GetXaxis()->SetTitle( " [TeV]" );
    interpolation->GetYaxis()->SetTitle( "n(#Lambda)" );
    interpolation->Draw( "AL" );
    ++nPad;
  }

  TH1 * dummy = new TH1D( "dummy", "", 11, 0.5, 11.5 );
  dummy->SetMaximum( 4.5 );
  dummy->SetMinimum( -7.5 );
  dummy->SetXTitle( "#Lambda [TeV]" );
  dummy->SetYTitle( "[ #mu^{MC}(#Lambda) - #mu^{fit}(#Lambda) ] / #sigma^{MC}" );
  vector<string> labels = assign::list_of("#infty")("8.00")("7.00")
      ("6.00")("5.00")("4.00")("3.00")("1.50")
      ("1.00")("0.75")("0.50");
  for( int bin = 1; bin <= dummy->GetNbinsX(); ++bin ) {
    dummy->GetXaxis()->SetBinLabel( bin, labels[bin-1].c_str() );
  }
  dummy->GetXaxis()->SetLabelSize( 0.1 );
  dummy->GetXaxis()->LabelsOption("v");
  dummy->GetXaxis()->SetTitleOffset(1.6);

  TLine line( 0.5, 0, 11.5, 0 );
  line.SetLineWidth( 2 );
  line.SetLineStyle( 2 );
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
            index.push_back( n - i ); // remember reverse ordering for alpha
            double graphY = yArr[i];
            double funcY = func->Eval( xArr[i] );
            offset.push_back( ( graphY - funcY ) / eyArr[i] );
          }
          TGraph * fitResult = new TGraph( n, &index[0], &offset[0] );
          _fitResults.push_back( fitResult );
          fitResult->SetMarkerColor( kRed );
          fitResult->SetTitle( title.c_str() );
          dummy->Draw("AXIS");
          line.Draw();
          fitResult->Draw( "P" );
          _latex.DrawLatex( 4, -5.0, ("#chi = " + CastToString<double>( chi )).c_str() );
          ++nPad;
        }

  typedef map< double, vector< double > > ParMap_t;
  ParMap_t params = pdf.pdfFitParams( mjj );
  vector< double > chis;
  vector< double > qcdPars;
  vector< double > qciPars;
  vector< double > interferencePars;
  foreach( ParMap_t::value_type& x, params )
  {
    chis.push_back( x.first );
    // see notes pg 126 for explanation
    double a0 = x.second[0];
    double a1 = x.second[1] / pow( 3., 2 );
    double a2 = x.second[2] / pow( 3., 4 );
    double chiSum = pdf.sumOverChi( pow( 3., -4 ) );
    double e2 = ( -( a1 - chiSum ) - sqrt( pow( a1 - chiSum, 2 ) - 4 * a0 * a2 ) ) / ( 2 * a2 );
    qcdPars.push_back( a0 / e2 );
    qciPars.push_back( a1 );
    interferencePars.push_back( a2 * e2 );
    //cout << "PDFMon chi = " << chis.back() << '\n' << "       qcd = " << qcdPars.back() << '\n' << "       qci = "
    //     << qciPars.back() << '\n' << "       int = " << interferencePars.back() << '\n';
  }

  {
    _parameterCanvas.cd();
    TGraph * qcd = new TGraph( chis.size(), &chis[0], &qcdPars[0] );
    qcd->SetLineColor( kBlack );
    qcd->SetLineWidth( 2 );
    qcd->SetName( "QCD" );
    _parameters.push_back( qcd );
    TGraph * qci = new TGraph( chis.size(), &chis[0], &qciPars[0] );
    qci->SetLineColor( kBlue );
    qci->SetLineWidth( 2 );
    qci->SetName( "CI" );
    _parameters.push_back( qci );
    TGraph * interference = new TGraph( chis.size(), &chis[0], &interferencePars[0] );
    interference->SetLineColor( kRed );
    interference->SetLineWidth( 2 );
    interference->SetName( "Interference" );
    _parameters.push_back( interference );
    TH2D* dummy = new TH2D( "paramMonDummy", "", 11, 0.9, 26, 100, -0.5, 1.5 );
    dummy->Draw();
    qcd->Draw( "LSAME" );
    qci->Draw( "LSAME" );
    interference->Draw( "LSAME" );
  }

  _pdfCanvas.Print( ( _folder + "PDFs" + _ext ).c_str() );
  _interpolCanvas.Print( ( _folder + "Interpolation" + _ext ).c_str() );
  _fitResultCanvas.Print( ( _folder + "FitResults" + _ext ).c_str() );
  _parameterCanvas.Print( ( _folder + "Parameters" + _ext ).c_str() );

}

