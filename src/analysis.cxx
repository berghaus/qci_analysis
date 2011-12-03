#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <stdexcept>
#include <limits>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/checked_delete.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/math/distributions/normal.hpp>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
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
using namespace boost::assign;
using namespace boost::math;

template< class T >
T quantile( const vector< T >&, const double& );

int main( int argc, char* argv[] ) {

  // TODO: give these as command line controls
  int nPE = 1000;
  double nBinsScale = 210;
  double minScale = 4.;
  double maxScale = 8.2;
  double deltaScale = ( maxScale - minScale ) / nBinsScale;
  vector< vector< double > > scales;
  for( double scale = minScale; scale < maxScale; scale += deltaScale )
    scales.push_back( vector< double >( 1, scale ) );

  // for GUI;
  TApplication theApp( "Analysis", &argc, argv );
  SetAtlasStyle();
  TGClient windowClient;
  const TGWindow * rootWindow = windowClient.GetRoot();
  ControlFrame * control = new ControlFrame( rootWindow, 350, 80 );

  TFile * pdfFile = TFile::Open( "~/docs/comp/analysis/kFactor.root", "READ" );
  TFile * dataFile = TFile::Open( "~/docs/comp/analysis/data.root", "READ" );
  TH1 * dataHist = (TH1*) dataFile->Get( "Chi_2000-to-7000all" );
  Experiment data( *dataHist );

  data.plot();
  cout << "integral(data) = " << data.integral() << endl;

  normal myGaus;
  // Set up PDF for data

  try {

    PDF * pdf = new PDF( pdfFile, data.integral() );
    Neg2LogLikelihoodRatio dataLikelihoodRatio( &data, pdf, 0. );

    foreach( const vector< double >& scale, scales )
      dataLikelihoodRatio( scale );

    TestStatMonitor tm( -1., "figures/", ".png" );
    for( int i = 0; i < 10; ++i ) {
      dataLikelihoodRatio.accept( tm );
      dataLikelihoodRatio.denominator().accept( tm );
    }
    tm.finalize();

    // Monitor data PDF
    PDFMonitor pdfMon;
    pdf->accept( pdfMon );

    PseudoExperimentFactory peFactory( pdf, data );

    typedef map< double, vector< PseudoExperiment > > peMap_t;

    // TODO: Encapsulate these in a class
    vector< double > scales;
    vector< double > observed;
    vector< double > expected;
    vector< double > e02;
    vector< double > e98;
    vector< double > e16;
    vector< double > e84;
    vector< double > ex;

    vector< PseudoExperiment > errorBandPEs = peFactory.build( 0., nPE );
    vector< Neg2LogLikelihoodRatio* > errorBandLRs;
    errorBandLRs.reserve( errorBandLRs.size() );
    foreach( const PseudoExperiment& pe, errorBandPEs )
    {
      PDF * pePDF = new PDF( pdf->pdfFitParams(), pe.integral() );
      errorBandLRs.push_back( new Neg2LogLikelihoodRatio( &pe, pePDF, 0. ) );
    }

    PValueTest pValueTest( 0., errorBandLRs );
    double pValue = pValueTest( dataLikelihoodRatio );
    cout << " * pvalue( Lambda = " << 0. << ") = " << pValue << endl;

    pValueTest.finalize();

    int scaleBin = 0;
    for( double scale = minScale; scale < maxScale; scale += deltaScale ) {
      double alpha = pow( scale, -4 );
      vector< PseudoExperiment > pes = peFactory.build( alpha, nPE );

      vector< Neg2LogLikelihoodRatio* > likelihoodRatios;
      likelihoodRatios.reserve( pes.size() );
      foreach( const PseudoExperiment& pe, pes )
      {
        PDF * pePDF = new PDF( pdf->pdfFitParams(), pe.integral() );
        likelihoodRatios.push_back( new Neg2LogLikelihoodRatio( &pe, pePDF, alpha ) );
      }

      PValueTest pv( alpha, likelihoodRatios ); // = *pValueTest;
      vector< double > par( 1, alpha );

      double pValue = pv( dataLikelihoodRatio );
      if ( !( scaleBin % 10 ) ) {
        pv.finalize();
      }
      ++scaleBin;

      vector< double > ebPValues;
      ebPValues.reserve( nPE );
      foreach( Neg2LogLikelihoodRatio* l, errorBandLRs )
        ebPValues.push_back( pv( *l ) );
      sort( ebPValues.begin(), ebPValues.end() );

      double sigmaN2 = quantile( ebPValues, cdf( myGaus, -2. ) );
      double sigmaN1 = quantile( ebPValues, cdf( myGaus, -1. ) );
      double exp = quantile( ebPValues, cdf( myGaus, 0. ) );
      double sigmaP1 = quantile( ebPValues, cdf( myGaus, 1. ) );
      double sigmaP2 = quantile( ebPValues, cdf( myGaus, 2. ) );

      scales.push_back( scale );
      observed.push_back( pValue );
      expected.push_back( exp );
      e02.push_back( fabs( exp - sigmaN2 ) );
      e98.push_back( fabs( exp - sigmaP2 ) );
      e16.push_back( fabs( exp - sigmaN1 ) );
      e84.push_back( fabs( exp - sigmaP1 ) );
      ex.push_back( 0. );

      // clean up likelihoods for this alpha
      for_each( likelihoodRatios.begin(), likelihoodRatios.end(), boost::checked_deleter< Neg2LogLikelihoodRatio >() );

    }
    // clean up error band Likelihoods
    for_each( errorBandLRs.begin(), errorBandLRs.end(), boost::checked_deleter< Neg2LogLikelihoodRatio >() );

    TGraph observedGraph( scales.size(), &scales[0], &observed[0] );
    observedGraph.SetLineWidth( 2 );
    observedGraph.SetLineColor( kRed );

    TGraphAsymmErrors expected1Sigma( nBinsScale, &scales[0], &expected[0], &ex[0], &ex[0], &e16[0], &e84[0] );
    expected1Sigma.SetFillStyle( 1001 );
    expected1Sigma.SetFillColor( kBlue );
    expected1Sigma.SetLineStyle( 2 );
    expected1Sigma.SetLineWidth( 2 );

    TGraphAsymmErrors expected2Sigma( nBinsScale, &scales[0], &expected[0], &ex[0], &ex[0], &e02[0], &e98[0] );
    expected2Sigma.SetFillStyle( 1001 );
    expected2Sigma.SetFillColor( kYellow );
    expected2Sigma.SetLineColor( kYellow );
    expected2Sigma.SetLineWidth( 2 );

    TMultiGraph exclusion;
    exclusion.Add( &expected2Sigma, "3" );
    exclusion.Add( &expected1Sigma, "3" );
    exclusion.Add( &expected1Sigma, "LX" );
    exclusion.Add( &observedGraph, "L" );

    TCanvas obsCanvas( "observedCanvas", "", 800, 800 );
    obsCanvas.cd();
    TH2D dummy( "dummy", "", 100, 2., 8., 100, 0., 1. );
    dummy.SetXTitle( "#Lambda [TeV]" );
    dummy.SetYTitle( "CL_{s+b}" );
    dummy.Draw();
    exclusion.Draw();
    dummy.Draw( "ASAME" );

    TCanvas obsCanvas2( "observedCanvas2", "", 800, 800 );
    obsCanvas2.cd();
    TH2D dummy2( "dummy2", "", 100, 2., 8, 100, 0.04, 0.1 );
    dummy2.SetXTitle( "#Lambda [TeV]" );
    dummy2.SetYTitle( "CL_{s+b}" );
    dummy2.Draw();
    exclusion.Draw();
    dummy2.Draw( "ASAME" );

    theApp.Run( kTRUE );

  } catch ( exception& e ) {
    // print exception to console
    cout << "caught exception:\n" << e.what() << endl;

    // give some time to look at problems
    theApp.Run( kTRUE );

  }

  pdfFile->Close();
  dataFile->Close();

  return 0;

}

// -----
template< class T >
T quantile( const vector< T >& v, const double& q ) {
  typename vector< T >::size_type index = v.size() * q;
  return v.at( index );
}

