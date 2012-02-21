/*
 * postProcess.cxx
 *
 *  Created on: 2012-01-11
 *      Author: frank
 */
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include <boost/checked_delete.hpp>
#include <boost/foreach.hpp>
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

#include <TClass.h>
#include <TCollection.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TKey.h>

#include "PostProcessCL.hpp"
#include "PValueTest.hpp"
#include "AtlasStyle.hpp"
#include "Effect.hpp"

#define foreach BOOST_FOREACH
#define ERROR_NO_SIGNAL_INPUT 1
#define ERROR_SIGNAL_BACKGROUND_MISMATCH 2
#define ERROR_NO_PDF 3
#define ERROR_NO_DATA 4

using namespace std;
namespace po = boost::program_options;

void ReadPValueTestsFromFile( const vector< string >&, vector< PValueTest >& );
vector< TDirectoryFile* > GetDirs( const TFile* );
vector< TH1* > GetHists( const TFile* );
template< class T > bool compByName( const T* x, const T* y ) {
  return string( x->GetName() ) < string( y->GetName() );
}

int main( int argc, char* argv[] ) {

  //---------------------------------------------------------------------------
  // process command line options
  // Declare the supported options.
  po::options_description desc( "Allowed options" );
  desc.add_options()( "help,h", "print this help message" )( "nPE,n", po::value< int >(),
                                                             "number of pseudo-experiments to run on each alpha" )(
      "sigInputFiles,s", po::value< vector< string > >()->multitoken(),
      "list of input files containing signal likelihood distributions" )(
      "bkgInputFiles,b", po::value< vector< string > >()->multitoken(),
      "list of input files containing background likelihood distributions" )( "outDir,o", po::value< string >(),
                                                                              "output directory for plots" )(
      "data,d", po::value< string >(), "ROOT file containing data event distribution" )(
      "pdf,p", po::value< string >(), "ROOT file containing expected event distributions" )(
          "stochastic", "include error due to limited statistics in QCD and QCI MC" )(
          "jes", po::value< string >(), "include error due to jet energy scale uncertainty described in given root file" )(
          "jer", po::value< string >(), "include error due to jet p_T resolution described in given root file" );

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

  vector< string > sigFileNames;
  if ( vm.count( "sigInputFiles" ) ) {
    sigFileNames = vm["sigInputFiles"].as< vector< string > >();
    cout << "reading signal likelihood distributions from:\n";
    foreach( const string& fn, sigFileNames )
      cout << "   * " << fn << "\n";
  } else {
    cout << "no signal input given. Aborting.\n";
    return ERROR_NO_SIGNAL_INPUT;
  }
  vector< PValueTest > sigLLDistributions;
  ReadPValueTestsFromFile( sigFileNames, sigLLDistributions );

  vector< string > bkgFileNames;
  bool doCLs = true;
  if ( vm.count( "bkgInputFiles" ) ) {
    bkgFileNames = vm["bkgInputFiles"].as< vector< string > >();
    cout << "reading background likelihood distributions from:\n";
    foreach( const string& fn, bkgFileNames )
      cout << "   * " << fn << "\n";
  } else {
    cout << "no background input given. disabeling CL_s computation.\n";
    doCLs = false;
  }
  vector< PValueTest > bkgLLDistributions;
  ReadPValueTestsFromFile( bkgFileNames, bkgLLDistributions );

  string folder = "./";
  if ( vm.count( "outDir" ) ) {
    folder = vm["outDir"].as< string >();
  }
  cout << "directing figures to: " << folder << "\n";

  string dataFileName;
  if ( vm.count( "data" ) ) {
    dataFileName = vm["data"].as< string >();
  } else {
    cout << "No data event distribution supplied. Aborting.\n";
    return ERROR_NO_DATA;
  }
  // read in data
  TFile * dataFile = TFile::Open( dataFileName.c_str(), "READ" );
  vector< TH1* > dataHists = GetHists( dataFile );
  TH1 * dataHist = dataHists.back();
  Experiment data( *dataHist );

  string pdfFileName;
  if ( vm.count( "pdf" ) ) {
    pdfFileName = vm["pdf"].as< string >();
  } else {
    cout << "No predicted event distributions supplied. Aborting.\n";
    return ERROR_NO_PDF;
  }
  // read in PDF
  TFile * pdfFile = TFile::Open( pdfFileName.c_str(), "READ" );
  vector< TDirectoryFile* > pdfDirs = GetDirs( pdfFile );
  TDirectoryFile* pdfDir = pdfDirs.back();
  Prediction * pdf = new Prediction( pdfDir, data.integral() );

  if ( vm.count( "stochastic" ) ) {
    cout << "including errors arising from limited statistics in QCD and QCI MC\n";
    Statitical_Effect * statEff = new Statitical_Effect( pdf->pdfFit( "PredictionFunctionForError" ),
                                                         pdf->covarianceMaticies() );
    Effect * eff = dynamic_cast< Effect* >( statEff );
    if ( eff ) pdf->addEffect( eff );
    else cout << "failed to downcast Statitical_Effect to Effect\n";

  }

  if ( vm.count( "jes" ) ) {
    cout << "including errors arising from jet energy scale uncertainty\n";
    string jesErrorFileName = vm["jes"].as< string >();
    JES_Systematic_Effect * sysEff = new JES_Systematic_Effect( jesErrorFileName );
    Effect * eff = dynamic_cast< Effect* >( sysEff );
    if ( eff ) pdf->addEffect( eff );
    else cout << "failed to downcast JES_Systematic_Effect to Effect\n";
  }


  if ( vm.count( "jer" ) ) {
    cout << "including errors arising from jet p_T Resolution\n";
    string jerErrorFileName = vm["jer"].as< string >();
    JER_Systematic_Effect * sysEff = new JER_Systematic_Effect( jerErrorFileName );
    Effect * eff = dynamic_cast< Effect* >( sysEff );
    if ( eff ) pdf->addEffect( eff );
    else cout << "failed to downcast JER_Systematic_Effect to Effect\n";
  }

  //---------------------------------------------------------------------------
  SetAtlasStyle();

  // sanity check for read in of likelihood distributions
  if ( doCLs && sigLLDistributions.size() != bkgLLDistributions.size() ) {
    cerr << "number of singal and background likelihood distributions does not match!\n";
    return ERROR_SIGNAL_BACKGROUND_MISMATCH;
  }
  // set up data likelihood distribution
  Neg2LogLikelihoodRatio dataLikelihoodRatio( &data, pdf, 0. );
  // --- make sure we get something reasonable across interesting scale values
  for( double scale = 2.; scale < 8.; scale += 0.1 )
    dataLikelihoodRatio( vector< double >( 1, scale ) );

  PseudoExperimentFactory peFactory( pdf, data );

  // set up likelihood ratios for error bands
  vector< PseudoExperiment > errorBandPEs = peFactory.build( 0., nPE );
  vector< Neg2LogLikelihoodRatio* > errorBandLRs;
  errorBandLRs.reserve( errorBandPEs.size() );
  foreach( const PseudoExperiment& pe, errorBandPEs )
  {
    Prediction * pePDF = new Prediction( *pdf );
    Neg2LogLikelihoodRatio * l = new Neg2LogLikelihoodRatio( &pe, pePDF, 0. );
    for( double scale = 2.; scale < 8.; scale += 0.1 )
      ( *l )( vector< double >( 1, scale ) );
    errorBandLRs.push_back( l );
  }

  // find p-value and limits
  if ( doCLs ) {
    PostProcessCL pp( sigLLDistributions, bkgLLDistributions, errorBandLRs, &dataLikelihoodRatio );
    pp.proc();
    pp.plot( folder );
    pp.print();
  } else {
    PostProcessCL pp( sigLLDistributions, errorBandLRs, &dataLikelihoodRatio );
    pp.proc();
    pp.plot( folder );
    pp.print();
  }

  // clean up error bar likelihoods
  for_each( errorBandLRs.begin(), errorBandLRs.end(), boost::checked_deleter< Neg2LogLikelihoodRatio >() );

  return 0;
}

void ReadPValueTestsFromFile( const vector< string >& names, vector< PValueTest >& pvs ) {

  foreach( const string& fn, names )
  {
    ifstream file;
    file.open( fn.c_str(), ios::binary );
    PValueTest buffer;
    while ( file >> buffer )
      pvs.push_back( buffer );
    file.close();
  }

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
