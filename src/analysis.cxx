#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <stdexcept>
#include <ctime>
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
#define ERROR_NO_PDF 3
#define ERROR_NO_DATA 4

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

  //---------------------------------------------------------------------------
  // process cmd opts
  // Declare the supported options.
  po::options_description desc( "Allowed options" );
  desc.add_options()( "help,h", "print this help message" )
      ( "nPE,n", po::value< int >(), "number of pseudo-experiments to run on each alpha" )
      ( "scales,s", po::value< vector< double > >()->multitoken(), "list of contact interaction scale values to run on (in TeV)" )
      ( "jobID,j", po::value< int >(), "PBS job ID for output naming" )
      ( "outDir,o", po::value< string >(), "output directory for likelihood disctributions" )
      ( "pdf,p", po::value< string >(), "ROOT file containing expected event distributions" )
      ( "data,d", po::value< string >(), "ROOT file containing data event distribution" );

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

  vector< double > scales;
  if ( vm.count( "scales" ) ) {
    scales = vm["scales"].as< vector< double > >();
    sort( scales.begin(), scales.end() );

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

  string outDir;
  if ( vm.count( "outDir" ) ) {
    outDir = vm["outDir"].as< string >();
    cout << "directing output to: " << outDir << "\n";
  }

  string pdfFileName;
  if ( vm.count( "pdf" ) ) {
    pdfFileName = vm["pdf"].as< string >();
  } else {
    cout << "No predicted event distributions supplied. Aborting.\n";
    return ERROR_NO_PDF;
  }

  string dataFileName;
  if ( vm.count( "data" ) ) {
    dataFileName = vm["data"].as< string >();
  } else {
    cout << "No data event distribution supplied. Aborting.\n";
    return ERROR_NO_DATA;
  }
  //---------------------------------------------------------------------------


  try {
    // Set up PDF to run

    // read in data file
    TFile * dataFile = TFile::Open( dataFileName.c_str(), "READ" );
    vector< TH1* > dataHists = GetHists( dataFile );
    TH1 * dataHist = dataHists.back();
    Experiment data( *dataHist );
    data.plot();

    // read in our PDF from file
    TFile * pdfFile = TFile::Open( pdfFileName.c_str(), "READ" );
    vector< TDirectoryFile* > pdfDirs = GetDirs( pdfFile );
    TDirectoryFile* pdfDir = pdfDirs.back();


    PDF * pdf = new PDF( pdfDir, data.integral() );
    Neg2LogLikelihoodRatio dataLikelihoodRatio( &data, pdf, 0. );

    // --- make sure we get something reasonable across interesting scale values
    for( double scale = 0.5; scale < 10.; scale += 0.1 )
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

    PseudoExperimentFactory peFactory( pdf, data, time( 0 ) );

    // create PEs for background likelihood distribution
    vector< PseudoExperiment > bgPEs = peFactory.build( 0., nPE );
    vector< Neg2LogLikelihoodRatio* > bgLikelihoodRatios;
    bgLikelihoodRatios.reserve( bgPEs.size() );
    foreach( const PseudoExperiment& pe, bgPEs )
    {
      PDF * pePDF = new PDF( *pdf );
      pePDF->nData( pe.integral() );
      Neg2LogLikelihoodRatio * l = new Neg2LogLikelihoodRatio( &pe, pePDF, 0. );
      for( double scale = 0.5; scale < 10.; scale += 0.1 )
            (*l)( vector< double >( 1, scale ) );
      bgLikelihoodRatios.push_back( l );
    }

    // --- open output files for likelihood distributions
    ofstream signalOutFile, bkgrndOutFile;
    signalOutFile.open( ( outDir + "/signalOutput-" + jobID + ".bin" ).c_str(), ios::binary );
    bkgrndOutFile.open( ( outDir + "/bkgrndOutput-" + jobID + ".bin" ).c_str(), ios::binary );

    int scaleBin = 0;
    foreach( double scale, scales )
    {
      double alpha = pow( scale, -4 );

      // create signal PEs
      vector< PseudoExperiment > sigPEs = peFactory.build( alpha, nPE );
      vector< Neg2LogLikelihoodRatio* > sigLikelihoodRatios;
      sigLikelihoodRatios.reserve( sigPEs.size() );
      foreach( const PseudoExperiment& pe, sigPEs )
      {
        PDF * pePDF = new PDF( *pdf );
        pePDF->nData( pe.integral() );
        Neg2LogLikelihoodRatio * l = new Neg2LogLikelihoodRatio( &pe, pePDF, alpha );
        for( double s = 0.5; s < 10.; s += 0.1 )
          ( *l )( vector< double >( 1, s ) );
        sigLikelihoodRatios.push_back( l );

      }

      // -------
      // CL_s+b
      PValueTest signalPlusBackgroundPValue( alpha, sigLikelihoodRatios );
      signalOutFile << signalPlusBackgroundPValue << endl;

      // -----------
      // CL_s
      PValueTest backgroundPValue( alpha, bgLikelihoodRatios ); // = *pValueTest;
      bkgrndOutFile << backgroundPValue << endl;

      // ----------
      // monitoring
      if ( !( scaleBin % 1000 ) ) {
        signalPlusBackgroundPValue.finalize();
        //backgroundPValue.finalize();
      }
      ++scaleBin;

    }

    // clean up background Likelihoods
    for_each( bgLikelihoodRatios.begin(), bgLikelihoodRatios.end(),
              boost::checked_deleter< Neg2LogLikelihoodRatio >() );

    // close output files
    signalOutFile.close();
    bkgrndOutFile.close();

    // close data and PDF input files
    pdfFile->Close();
    dataFile->Close();

  } catch ( exception& e ) {
    // print exception to console
    cout << "caught exception:\n" << e.what() << endl;

  }

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

