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
#include <TGraphErrors.h>

#include "Likelihood.hpp"
#include "Prediction.hpp"
#include "PseudoExperiment.hpp"
#include "PValueTest.hpp"
#include "AtlasStyle.hpp"
#include "ControlFrame.hpp"
#include "CertaintyLevel.hpp"
#include "Effect.hpp"

#include "PredictionMonitor.hpp"
#include "TestStatMonitor.hpp"

#define ERROR_NO_SCALE_VALUE 1
#define ERROR_NO_PDF 3
#define ERROR_NO_DATA 4

using boost::lexical_cast;
using namespace std;
using namespace boost::assign;
namespace po = boost::program_options;

map< double, TDirectoryFile* > GetDirs( const TFile* );
map< double, TH1* > GetHists( const TFile* );
template< class T > bool compByName( const T* x, const T* y ) {
  return string( x->GetName() ) < string( y->GetName() );
}

void PlotPredictionError( Experiment&, Prediction *, Statitical_Effect * );

int main( int argc, char* argv[] ) {

  //---------------------------------------------------------------------------
  // process cmd opts
  // Declare the supported options.
  po::options_description desc( "Allowed options" );
  desc.add_options()( "help,h", "print this help message" )( "nPE,n", po::value< int >(),
                                                             "number of pseudo-experiments to run on each alpha" )(
      "scales,s", po::value< vector< double > >()->multitoken(), "list of contact interaction scale values (in TeV)" )(
      "mjjs,m", po::value< vector< double > >()->multitoken(), "mjj bins to consider (use mjj max)" )(
      "jobID,j", po::value< int >(), "PBS job ID for output naming" )(
      "outDir,o", po::value< string >(), "output directory for likelihood disctributions" )(
      "prediction,p", po::value< string >(), "ROOT file containing expected event distributions" )(
      "data,d", po::value< string >(), "ROOT file containing data event distribution" )(
      "stochastic", "include error due to limited statistics in QCD and QCI MC" )(
      "jes", po::value< string >(), "include error due to jet energy scale uncertainty described in given root file" )(
      "jer", po::value< string >(), "include error due to jet p_T resolution described in given root file" )(
      "pdf", po::value< string >(), "include error due to PDF sets as described in given root file" )(
      "appl", po::value< vector< string > >()->multitoken(), "include error due to scale choices as described in given grid files" )(
      "figures,f", po::value< string >(), "directory for output figures" );

  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, desc ), vm );
  po::notify( vm );

  if ( vm.count( "help" ) ) {
    cout << desc << "\n";
    return 0;
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

  vector< double > mjjs( 1, 7000. );
  if ( vm.count( "mjjs" ) ) {
    mjjs = vm["mjjs"].as< vector< double > >();
    sort( mjjs.begin(), mjjs.end() );

    cout << "Running over " << mjjs.size() << " scale values:\n";
    foreach( double mjj, mjjs )
      cout << "   - " << mjj << '\n';
  } else {
    cout << "No scale to run on given, defaulting to top mjj bin." << endl;
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

  string dataFileName;
  if ( vm.count( "data" ) ) {
    dataFileName = vm["data"].as< string >();
  } else {
    cout << "No data event distribution supplied. Aborting.\n";
    return ERROR_NO_DATA;
  }
  // read in data file
  TFile * dataFile = TFile::Open( dataFileName.c_str(), "READ" );
  typedef map< double, TH1* > th1Map_t;
  th1Map_t dataHists = GetHists( dataFile );
  th1Map_t::iterator dataIt = dataHists.begin();
  while ( dataIt != dataHists.end() ) {
    if ( find( mjjs.begin(), mjjs.end(), dataIt->first ) == mjjs.end() ) dataHists.erase( dataIt++ );
    else ++dataIt;
  }
  Experiment data( dataHists );
  //data.plot();

  string predictionFileName;
  if ( vm.count( "prediction" ) ) {
    predictionFileName = vm["prediction"].as< string >();
  } else {
    cout << "No predicted event distributions supplied. Aborting.\n";
    return ERROR_NO_PDF;
  }
  // read in our PDF from file
  TFile * pdfFile = TFile::Open( predictionFileName.c_str(), "READ" );
  typedef map< double, TDirectoryFile* > tDirFileMap_t;
  tDirFileMap_t pdfDirs = GetDirs( pdfFile );
  tDirFileMap_t::iterator predIt = pdfDirs.begin();
  while ( predIt != pdfDirs.end() ) {
    if ( find( mjjs.begin(), mjjs.end(), predIt->first ) == mjjs.end() ) pdfDirs.erase( predIt++ );
    else ++predIt;
  }
  // Set up PDF to run
  Prediction * pdf = new Prediction( pdfDirs, data );

  if ( vm.count( "stochastic" ) ) {
    cout << "including errors arising from limited statistics in QCD and QCI MC\n";
    Statitical_Effect * statEff = new Statitical_Effect( pdf->pdfFit( "PredictionFunctionForError" ),
                                                         pdf->covarianceMaticies() );
    Effect * eff = dynamic_cast< Effect* >( statEff );
    if ( eff ) pdf->addEffect( eff );
    else cout << "failed to downcast Statitical_Effect to Effect\n";

    PlotPredictionError( data, pdf, statEff );

  }

  if ( vm.count( "jes" ) ) {
    cout << "including errors arising from jet energy scale uncertainty\n";
    string jesErrorFileName = vm["jes"].as< string >();
    JES_Effect * sysEff = new JES_Effect( jesErrorFileName );
    Effect * eff = dynamic_cast< Effect* >( sysEff );
    if ( eff ) pdf->addEffect( eff );
    else cout << "failed to downcast JES_Effect to Effect\n";
  }

  if ( vm.count( "jer" ) ) {
    cout << "including errors arising from jet p_T Resolution\n";
    string jerErrorFileName = vm["jer"].as< string >();
    JER_Effect * sysEff = new JER_Effect( jerErrorFileName );
    Effect * eff = dynamic_cast< Effect* >( sysEff );
    if ( eff ) pdf->addEffect( eff );
    else cout << "failed to downcast JER_Effect to Effect\n";
  }

  if( vm.count( "pdf" ) ) {
    cout << "including errors arising from PDF fit\n";
    string pdfErrorFileName = vm["pdf"].as< string > ();
    PDF_Effect * sysEff = new PDF_Effect( pdfErrorFileName );
    Effect * eff = dynamic_cast< Effect* > ( sysEff );
    if( eff )
      pdf->addEffect( eff );
    else
      cout << "failed to downcast PDF_Effect to Effect\n";
  }

  if( vm.count( "appl" ) ) {
    cout << "including errors arising from scale choice\n";
    vector< string > gridNames = vm["appl"].as< vector< string > > ();
    Scale_Effect * sysEff = new Scale_Effect( gridNames );
    Effect * eff = dynamic_cast< Effect* > ( sysEff );
    if( eff )
      pdf->addEffect( eff );
    else
      cout << "failed to downcast Scale_Effect to Effect\n";
  }


  string figureDir = "./";
  if ( vm.count( "figures" ) ) {
    figureDir = vm["figures"].as< string >();
  }
  cout << "directing figures to " << figureDir << "\n";
  //---------------------------------------------------------------------------
  SetAtlasStyle();

  try {

    Neg2LogMaximumLikelihoodRatio dataLikelihoodRatio( &data, pdf, 0. );

    // --- make sure we get something reasonable across interesting scale values
    for( double scale = 0.5; scale < 10.; scale += 0.1 )
      dataLikelihoodRatio( vector< double >( 1, scale ) );

    TestStatMonitor tm( -1., figureDir + "/Likelihood/", ".pdf" );
    for( int i = 0; i < 10; ++i ) {
      dataLikelihoodRatio.accept( tm );
      dataLikelihoodRatio.denominator().accept( tm );
    }
    tm.finalize();

    // Monitor data PDF
    PredictionMonitor pdfMon( figureDir + "/PDF/", ".pdf" );
    pdf->accept( pdfMon );

    PseudoExperimentFactory peFactory( pdf, data );

    // create PEs for background likelihood distribution
    vector< PseudoExperiment > bgPEs = peFactory.build( 0., nPE );
    vector< Neg2LogMaximumLikelihoodRatio* > bgLikelihoodRatios;
    bgLikelihoodRatios.reserve( bgPEs.size() );
    foreach( const PseudoExperiment& pe, bgPEs )
    {
      Prediction * pePDF = new Prediction( *pdf );
      pePDF->nData( dynamic_cast< const Experiment& >( pe ) );
      Neg2LogMaximumLikelihoodRatio * l = new Neg2LogMaximumLikelihoodRatio( &pe, pePDF, 0. );
      for( double scale = 0.5; scale < 10.; scale += 0.1 )
        ( *l )( vector< double >( 1, scale ) );
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
      vector< Neg2LogMaximumLikelihoodRatio* > sigLikelihoodRatios;
      sigLikelihoodRatios.reserve( sigPEs.size() );
      foreach( const PseudoExperiment& pe, sigPEs )
      {
        Prediction * pePDF = new Prediction( *pdf );
        pePDF->nData( dynamic_cast< const Experiment& >( pe ) );
        Neg2LogMaximumLikelihoodRatio * l = new Neg2LogMaximumLikelihoodRatio( &pe, pePDF, alpha );
        for( double s = 0.5; s < 10.; s += 0.1 )
          ( *l )( vector< double >( 1, s ) );
        sigLikelihoodRatios.push_back( l );

      }

      // -------
      // CL_s+b
      PValueTest<Neg2LogMaximumLikelihoodRatio> signalPlusBackgroundPValue( alpha, sigLikelihoodRatios );
      signalOutFile << signalPlusBackgroundPValue << endl;

      // -----------
      // CL_s
      PValueTest<Neg2LogMaximumLikelihoodRatio> backgroundPValue( alpha, bgLikelihoodRatios ); // = *pValueTest;
      bkgrndOutFile << backgroundPValue << endl;

      // ----------
      // monitoring
      if ( !( scaleBin % 1000 ) ) {
        signalPlusBackgroundPValue.finalize( figureDir + "/Likelihood/signal/" );
        backgroundPValue.finalize( figureDir + "/Likelihood/bkgrnd/" );
      }
      ++scaleBin;

    }

    // clean up background Likelihoods
    for_each( bgLikelihoodRatios.begin(), bgLikelihoodRatios.end(),
              boost::checked_deleter< Neg2LogMaximumLikelihoodRatio >() );

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

void PlotPredictionError( Experiment& d, Prediction * p, Statitical_Effect * e ) {

  SetAtlasStyle();

  p->setMjj( 7000. );

  int nBins = 1000;
  double chiMin = 2.;
  double chiMax = 8.;
  vector<double> scales;
  vector<double> y;
  vector<double> ex;
  vector<double> ey;

  TCanvas canvas( "fitErrorCan", "", 850, 1100 );
  canvas.Divide( 3, 4, 0, 0 );

  for ( int i = 0; i <= nBins; ++i ) {
    scales.push_back( chiMin + double(i) * ( chiMax - chiMin ) / double( nBins ) );
  }

  TH1D dummy("dummy","",2, chiMin, chiMax + 0.5 );
  dummy.SetMinimum(  1. );
  dummy.SetMaximum( 85. );
  dummy.SetXTitle( "#Lambda [TeV]" );
  dummy.SetYTitle( "#mu_{j}(#Lambda)" );

  for( int bin = 0; bin < d[7000.].chi().size(); ++bin ) {
    y.clear();
    ex.clear();
    ey.clear();

    double chi = Prediction::labelChi( d[7000.].chi( bin ) );
    foreach( double l, scales ) {
      double alpha = pow( l, -4 );
      double expectedN = p->interpolate( chi, alpha);
      y.push_back( expectedN );
      ex.push_back( ( chiMax - chiMin ) / double( nBins ) / 2. );
      ey.push_back( e->error(expectedN, l, chi, 7000.) );
    }

    canvas.cd( bin+1 );
    dummy.Draw("AXIS");
    TGraphErrors * g = new TGraphErrors( nBins, &scales[0], &y[0], &ex[0], &ey[0]);
    g->SetLineColor( kBlue );
    g->SetFillColor( kBlue );
    g->Draw("3");
    dummy.Draw("SAMEAXIS");
  }

  canvas.Print( "fitErrors.pdf" );

}
