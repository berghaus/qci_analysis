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
#include "CertaintyLevel.hpp"

#include "PDFMonitor.hpp"
#include "TestStatMonitor.hpp"

using namespace std;
using namespace boost::assign;


int main( int argc, char* argv[] ) {

  // TODO: give these as command line controls
  int nPE = 2000;
  double nBinsScale = 60;
  double minScale = 2.;
  double maxScale = 8.;
  double deltaScale = ( maxScale - minScale ) / nBinsScale;
  vector< vector< double > > scales;
  for( double scale = minScale; scale < maxScale; scale += deltaScale )
    scales.push_back( vector< double >( 1, scale ) );

  CertaintyLevel CL_sb("CL_{s+b}", nBinsScale, minScale, maxScale );
  CertaintyLevel CL_s("CL_{s}", nBinsScale, minScale, maxScale );

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

    pValueTest.finalize();

    
    vector< PseudoExperiment > bgPEs = peFactory.build( 0., nPE );
    vector< Neg2LogLikelihoodRatio* > bgLikelihoodRatios;
    bgLikelihoodRatios.reserve( bgPEs.size() );
    foreach( const PseudoExperiment& pe, bgPEs )
    {
      PDF * pePDF = new PDF( pdf->pdfFitParams(), pe.integral() );
      bgLikelihoodRatios.push_back( new Neg2LogLikelihoodRatio( &pe, pePDF, 0. ) );
    }    
    
    int scaleBin = 0;
    for( double scale = minScale; scaleBin < nBinsScale; scale += deltaScale ) {
      double alpha = pow( scale, -4 );
      vector< PseudoExperiment > pes = peFactory.build( alpha, nPE );

      vector< Neg2LogLikelihoodRatio* > likelihoodRatios;
      likelihoodRatios.reserve( pes.size() );
      foreach( const PseudoExperiment& pe, pes )
      {
        PDF * pePDF = new PDF( pdf->pdfFitParams(), pe.integral() );
        likelihoodRatios.push_back( new Neg2LogLikelihoodRatio( &pe, pePDF, alpha ) );
      }

      // -------
      // CL_s+b
      PValueTest signalPlusBackgroundPValue( alpha, likelihoodRatios ); // = *pValueTest;
      vector< double > par( 1, alpha );
      double clsb_observed = signalPlusBackgroundPValue( dataLikelihoodRatio );

      vector< double > clsb_expected;
      clsb_expected.reserve( nPE );
      foreach( Neg2LogLikelihoodRatio* l, errorBandLRs )
        clsb_expected.push_back( signalPlusBackgroundPValue( *l ) );
      sort( clsb_expected.begin(), clsb_expected.end() );

      CL_sb.add( scale, clsb_observed, clsb_expected );


      // -----------
      // CL_s
      PValueTest backgroundPValue( alpha, bgLikelihoodRatios ); // = *pValueTest;
      double cls_observed = signalPlusBackgroundPValue( dataLikelihoodRatio ) / backgroundPValue( dataLikelihoodRatio );

      vector< double > cls_expected;
      cls_expected.reserve( nPE );
      foreach( Neg2LogLikelihoodRatio* l, errorBandLRs )
        cls_expected.push_back( signalPlusBackgroundPValue( *l ) / backgroundPValue( *l ) );
      sort( cls_expected.begin(), cls_expected.end() );

      CL_s.add( scale, cls_observed, cls_expected );

      // ----------
      // monitoring
      if ( !( scaleBin % 10 ) ) {
        signalPlusBackgroundPValue.finalize();
        backgroundPValue.finalize();
      }
      ++scaleBin;

      // clean up likelihoods for this alpha
      for_each( likelihoodRatios.begin(), likelihoodRatios.end(), boost::checked_deleter< Neg2LogLikelihoodRatio >() );

    }
    // clean up error band Likelihoods
    for_each( errorBandLRs.begin(), errorBandLRs.end(), boost::checked_deleter< Neg2LogLikelihoodRatio >() );
    for_each( bgLikelihoodRatios.begin(), bgLikelihoodRatios.end(), boost::checked_deleter< Neg2LogLikelihoodRatio >() );

    CL_sb.plot();
    CL_s.plot();


    cout << " * pvalue( Lambda = " << 0. << ") = " << pValue << endl;
    cout << CL_sb << endl;
    cout << CL_s << endl;

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

