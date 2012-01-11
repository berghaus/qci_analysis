#include <algorithm>
#include <fstream>
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
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/config.hpp>
#include <boost/program_options/environment_iterator.hpp>
#include <boost/program_options/eof_iterator.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/version.hpp>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH2.h>
#include <TFile.h>
#include <TKey.h>
#include <TDirectoryFile.h>
#include <TClass.h>

#include "Likelihood.hpp"
#include "PDF.hpp"
#include "PseudoExperiment.hpp"
#include "PValueTest.hpp"
#include "AtlasStyle.hpp"
#include "ControlFrame.hpp"
#include "CertaintyLevel.hpp"

#include "PDFMonitor.hpp"
#include "TestStatMonitor.hpp"

#define ERROR_NO_SCALE_VALUE 1

using boost::lexical_cast;
using namespace std;
using namespace boost::assign;
namespace po = boost::program_options;

vector< TDirectoryFile* > GetDirs( const TFile* );
vector< TH1* > GetHists( const TFile* );
template< class T > bool compByName( const T* x, const T* y ) {
  return string( x->GetName() ) < string( y->GetName() );
}

int main( int argc, char* argv[] ) {

  // process cmd opts
  // Declare the supported options.
  po::options_description desc( "Allowed options" );
  desc.add_options()( "help,h", "print this help message" )( "nPE,n", po::value< int >(),
                                                             "number of pseudo-experiments to run on each alpha" )(
      "scales,s", po::value< vector< double > >()->multitoken(),
      "list of contact interaction scale values to run on (in TeV)" )( "jobID,j", po::value< int >(),
                                                                       "PBS job ID for output naming" );
  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, desc ), vm );
  po::notify( vm );

  if ( vm.count( "help" ) ) {
    cout << desc << "\n";
    return 1;
  }

  int nPE = 1000;
  if ( vm.count( "nPE" ) ) {
    cout << "Running " << vm["nPE"].as< int >() << " pseudo-experiments for each alpha.\n";
    nPE = vm["nPE"].as< int >();
  } else {
    cout << "Defaulting to using " << nPE << " pseudo-experiments for each alpha.\n";
  }

  // TODO: give these as command line controls
  double nBinsScale = 0.;
  double minScale = 0.;
  double maxScale = 0.;
  vector< double > scales;
  if ( vm.count( "scales" ) ) {
    scales = vm["scales"].as< vector< double > >();
    sort( scales.begin(), scales.end() );
    minScale = scales.front();
    maxScale = scales.back();
    nBinsScale = scales.size();
    cout << "Running over " << scales.size() << " scale values:\n";
    foreach( double scale, scales )
      cout << "   - " << scale << '\n';
  } else {
    cout << "No scale to run on given ... aborting." << endl;
    return ERROR_NO_SCALE_VALUE;
  }

  string jobID = "";
  if ( vm.count( "jobID" ) ) {
    int j = vm["jobID"].as< int >();
    cout << "PBS job ID: " << j << "\n";
    jobID = lexical_cast< string >( j );
  }

  CertaintyLevel CL_sb( "CL_{s+b}", nBinsScale, minScale, maxScale );
  CertaintyLevel CL_s( "CL_{s}", nBinsScale, minScale, maxScale );

  // for GUI;
//  TApplication theApp( "Analysis", &argc, argv );
//  SetAtlasStyle();
//  TGClient windowClient;
//  const TGWindow * rootWindow = windowClient.GetRoot();
//  ControlFrame * control = new ControlFrame( rootWindow, 350, 80 );

  // read in our files
  TFile * pdfFile = TFile::Open( "~/docs/comp/analysis/kFactor.root", "READ" );
  vector< TDirectoryFile* > pdfDirs = GetDirs( pdfFile );
  TDirectoryFile* pdfDir = pdfDirs.back();

  TFile * dataFile = TFile::Open( "~/docs/comp/analysis/data.root", "READ" );
  vector< TH1* > dataHists = GetHists( dataFile );
  TH1 * dataHist = dataHists.back();
  Experiment data( *dataHist );

  data.plot();

  // Set up PDF for data

  try {

    PDF * pdf = new PDF( pdfDir, data.integral() );
    Neg2LogLikelihoodRatio dataLikelihoodRatio( &data, pdf, 0. );

    foreach( const double& scale, scales )
      dataLikelihoodRatio( vector< double >( 1, scale ) );

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

    // -- file to ouput PValue test into
    ofstream signalOutFile, bkgrndOutFile;
    signalOutFile.open( ("signalOutput-" + jobID + ".bin").c_str(), ios::binary );
    bkgrndOutFile.open( ("bkgrndOutput-" + jobID + ".bin").c_str(), ios::binary );
    int scaleBin = 0;
    foreach( double scale, scales )
    {
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

      signalOutFile << signalPlusBackgroundPValue;

      // -----------
      // CL_s
      PValueTest backgroundPValue( alpha, bgLikelihoodRatios ); // = *pValueTest;
      bkgrndOutFile << backgroundPValue;
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
    for_each( bgLikelihoodRatios.begin(), bgLikelihoodRatios.end(),
              boost::checked_deleter< Neg2LogLikelihoodRatio >() );

    CL_sb.plot();
    CL_s.plot();

    cout << " * pvalue( Lambda = " << 0. << ") = " << pValue << endl;
    cout << CL_sb << endl;
    cout << CL_s << endl;

//    theApp.Run( kTRUE );
    signalOutFile.close();
    bkgrndOutFile.close();

  } catch ( exception& e ) {
    // print exception to console
    cout << "caught exception:\n" << e.what() << endl;

    // give some time to look at problems
//    theApp.Run( kTRUE );

  }

  pdfFile->Close();
  dataFile->Close();

  return 0;

}

vector< TDirectoryFile* > GetDirs( const TFile* file ) {

  vector< TDirectoryFile* > result;

  TIter nextKey( file->GetListOfKeys() );
  TKey * key = (TKey*) nextKey();

  do {

    if ( !key ) break;

    string name = key->GetName();
    if ( name.find( "0" ) == 0 ) continue;
    TObject * obj = key->ReadObj();

    if ( obj && obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      cout << "found: " << name << endl;
      result.push_back( (TDirectoryFile*) obj );
    }

  } while ( key = (TKey*) nextKey() );

  return result;

}

vector< TH1* > GetHists( const TFile* file ) {

  vector< TH1* > result;

  TIter nextKey( file->GetListOfKeys() );
  TKey * key = (TKey*) nextKey();

  do {

    if ( !key ) break;

    string name = key->GetName();
    if ( name.find( "Chi_" ) != 0 ) continue;
    TObject * obj = key->ReadObj();

    if ( obj && obj->IsA()->InheritsFrom( "TH1" ) ) {
      cout << "found: " << name << endl;
      result.push_back( (TH1*) obj );
    }

  } while ( key = (TKey*) nextKey() );

  return result;

}

