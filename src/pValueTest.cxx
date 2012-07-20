/*
 * pValueTest.cxx
 *
 *  Created on: 2012-07-19
 *      Author: frank
 */
#include <map>
#include <string>
#include <vector>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
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

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TKey.h>
#include <TDirectoryFile.h>
#include <TClass.h>
#include <TGraphErrors.h>

#include "Likelihood.hpp"
#include "Prediction.hpp"
#include "PseudoExperiment.hpp"
#include "PValueTest.hpp"
#include "AtlasStyle.hpp"


#define ERROR_NO_PDF 3
#define ERROR_NO_DATA 4

using namespace std;
using boost::lexical_cast;
namespace po = boost::program_options;

map< double, TDirectoryFile* > GetDirs( const TFile* );
map< double, TH1* > GetHists( const TFile* );
template< class T > bool compByName( const T* x, const T* y ) {
  return string( x->GetName() ) < string( y->GetName() );
}


int main( int argc, char* argv[] ) {

  // variables will be set by cmd line options
  int nPE = 1000;
  vector< double > mjjs( 1, 7000. );
  string dataFileName;
  string predictionFileName;
  string figureDir = "./";

  //---------------------------------------------------------------------------
  // process cmd opts
  // Declare the supported options.
  po::options_description desc( "Allowed options" );
  desc.add_options()( "help,h", "print this help message" )
      ( "nPE,n", po::value< int >(), "number of pseudo-experiments to run on each alpha" )
      ( "mjjs,m", po::value< vector< double > >()->multitoken(), "mjj bins to consider (use mjj max)" )
      ( "prediction,p", po::value< string >(), "ROOT file containing expected event distributions" )
      ( "data,d", po::value< string >(), "ROOT file containing data event distribution" )
      ( "figures,f", po::value< string >(), "directory for output figures" );

  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, desc ), vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    cout << desc << "\n";
    return 0;
  }

  if( vm.count( "nPE" ) ) {
    cout << "Running " << vm["nPE"].as< int > () << " pseudo-experiments for each alpha.\n";
    nPE = vm["nPE"].as< int > ();
  } else {
    cout << "Defaulting to using " << nPE << " pseudo-experiments for each alpha.\n";
  }

  if( vm.count( "mjjs" ) ) {
    mjjs = vm["mjjs"].as< vector< double > > ();
    sort( mjjs.begin(), mjjs.end() );

    cout << "Running over " << mjjs.size() << " scale values:\n";
    foreach ( double mjj, mjjs )
            cout << "   - " << mjj << '\n';
  } else {
    cout << "No scale to run on given, defaulting to top mjj bin." << endl;
  }

  if( vm.count( "data" ) ) {
    dataFileName = vm["data"].as< string > ();
  } else {
    cout << "No data event distribution supplied. Aborting.\n";
    return ERROR_NO_DATA;
  }
  // read in data file
  TFile * dataFile = TFile::Open( dataFileName.c_str(), "READ" );
  typedef map< double, TH1* > th1Map_t;
  th1Map_t dataHists = GetHists( dataFile );
  th1Map_t::iterator dataIt = dataHists.begin();
  while( dataIt != dataHists.end() ) {
    if( find( mjjs.begin(), mjjs.end(), dataIt->first ) == mjjs.end() )
      dataHists.erase( dataIt++ );
    else
      ++dataIt;
  }
  Experiment data( dataHists );
  //data.plot();

  if( vm.count( "prediction" ) ) {
    predictionFileName = vm["prediction"].as< string > ();
  } else {
    cout << "No predicted event distributions supplied. Aborting.\n";
    return ERROR_NO_PDF;
  }
  // read in our PDF from file
  TFile * pdfFile = TFile::Open( predictionFileName.c_str(), "READ" );
  typedef map< double, TDirectoryFile* > tDirFileMap_t;
  tDirFileMap_t pdfDirs = GetDirs( pdfFile );
  tDirFileMap_t::iterator predIt = pdfDirs.begin();
  while( predIt != pdfDirs.end() ) {
    if( find( mjjs.begin(), mjjs.end(), predIt->first ) == mjjs.end() )
      pdfDirs.erase( predIt++ );
    else
      ++predIt;
  }
  // Set up PDF to run
  Prediction * prediction = new Prediction( pdfDirs, data );

  if( vm.count( "figures" ) ) {
    figureDir = vm["figures"].as< string > ();
  }
  cout << "directing figures to " << figureDir << "\n";
  //---------------------------------------------------------------------------
  SetAtlasStyle();

  try {

      Neg2LogMaximumLikelihoodRatio dataLLR( &data, prediction, 0. );
      Neg2LogLikelihood_FCN dataLL( &data, prediction, 0. );

      // --- make sure we get something reasonable across interesting scale values
      // for( double scale = 0.5; scale < 10.; scale += 0.1 )
      //   dataLLR( vector< double >( 1, scale ) );

      PseudoExperimentFactory peFactory( prediction, data );

      // create PEs for background likelihood distribution
      vector< PseudoExperiment > pseudoExperiments = peFactory.build( 0., nPE );
      vector< Neg2LogMaximumLikelihoodRatio* > llrs;
      vector< Neg2LogLikelihood_FCN* > lls;
      llrs.reserve( pseudoExperiments.size() );
      lls.reserve( pseudoExperiments.size() );
      foreach( const PseudoExperiment& pe, pseudoExperiments )
      {
        Prediction * pePrediction = new Prediction( *prediction );
        pePrediction->nData( dynamic_cast< const Experiment& >( pe ) );

        Neg2LogLikelihood_FCN * peLL = new Neg2LogLikelihood_FCN( &pe, prediction, 0. );
        lls.push_back( peLL );

        // Neg2LogMaximumLikelihoodRatio * l = new Neg2LogMaximumLikelihoodRatio( &pe, pePrediction, 0. );
//        for( double scale = 0.5; scale < 10.; scale += 0.1 )
//          ( *l )( vector< double >( 1, scale ) );
//        llrs.push_back( l );


      }

      // -------
      // p-value from Poisson likelihood
      PValueTest<Neg2LogLikelihood_FCN> LLpValue( 0., lls );
      cout << "\n\n   * p-value = " <<  LLpValue( dataLL ) << "\n\n";
      LLpValue.finalize( figureDir );

      // -------
      // p-value from likelihood ratio
//      PValueTest<Neg2LogMaximumLikelihoodRatio> LLRpValue( 0., llrs );
//      cout << "\n\n   * p-value = " <<  LLRpValue( dataLLR ) << "\n\n";
//      LLRpValue.finalize( figureDir );


  } catch( ... ) {

    cout << "stumbled across an exception!\n";

  }

}

map< double, TDirectoryFile* > GetDirs( const TFile* file ) {

  map< double, TDirectoryFile* > result;

  TIter nextKey( file->GetListOfKeys() );
  TKey * key = (TKey*) nextKey();

  do {

    if ( !key ) break;

    string name = key->GetName();
    if ( name.find( "0" ) == 0 ) continue;
    TObject * obj = key->ReadObj();

    if ( obj && obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      cout << "found: " << name << endl;
      double mjjMax = lexical_cast< double >(
          name.substr( name.find( "-mjj-" ) + 5, name.size() - name.find( "-mjj-" ) - 8 ) );
      result.insert( make_pair( mjjMax, (TDirectoryFile*) obj ) );
    }

  } while ( key = (TKey*) nextKey() );

  return result;

}

map< double, TH1* > GetHists( const TFile* file ) {

  map< double, TH1* > result;

  TIter nextKey( file->GetListOfKeys() );
  TKey * key = (TKey*) nextKey();

  do {

    if ( !key ) break;

    string name = key->GetName();
    if ( name.find( "Chi_" ) != 0 ) continue;
    TObject * obj = key->ReadObj();

    if ( obj && obj->IsA()->InheritsFrom( "TH1" ) ) {
      cout << "found: " << name << endl;
      double mjjMax = lexical_cast< double >(
          name.substr( name.find( "-to-" ) + 4, name.size() - name.find( "-to-" ) - 7 ) );
      result.insert( make_pair( mjjMax, (TH1*) obj ) );
    }

  } while ( key = (TKey*) nextKey() );

  return result;

}

