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
#include <boost/format.hpp>
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

#include "Prediction.hpp"
#include "PseudoExperiment.hpp"
#include "AtlasStyle.hpp"
#include "Effect.hpp"

#include "PredictionMonitor.hpp"

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
void plotErrors( const Experiment&, const Prediction&, const Statitical_Effect&, const double& );

string g_outDir = "./";

int main( int argc, char* argv[] ) {

  //---------------------------------------------------------------------------
  // process cmd opts
  // Declare the supported options.
  po::options_description desc( "Allowed options" );
  desc.add_options()( "help,h", "print this help message" )(
      "scales,s", po::value< vector< double > >()->multitoken(), "list of contact interaction scale values (in TeV)" )(
      "mjjs,m", po::value< vector< double > >()->multitoken(), "mjj bins to consider (use mjj max)" )(
      "outDir,o", po::value< string >(), "output directory for likelihood disctributions" )(
      "pdf,p", po::value< string >(), "ROOT file containing expected event distributions" )(
      "data,d", po::value< string >(), "ROOT file containing data event distribution" );

  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, desc ), vm );
  po::notify( vm );

  if ( vm.count( "help" ) ) {
    cout << desc << "\n";
    return 0;
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

  if ( vm.count( "outDir" ) ) {
    g_outDir = vm["outDir"].as< string >();
    cout << "directing output to: " << g_outDir << "\n";
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

  string pdfFileName;
  if ( vm.count( "pdf" ) ) {
    pdfFileName = vm["pdf"].as< string >();
  } else {
    cout << "No predicted event distributions supplied. Aborting.\n";
    return ERROR_NO_PDF;
  }
  // read in our PDF from file
  TFile * pdfFile = TFile::Open( pdfFileName.c_str(), "READ" );
  typedef map< double, TDirectoryFile* > tDirFileMap_t;
  tDirFileMap_t pdfDirs = GetDirs( pdfFile );
  tDirFileMap_t::iterator predIt = pdfDirs.begin();
  while ( predIt != pdfDirs.end() ) {
    if ( find( mjjs.begin(), mjjs.end(), predIt->first ) == mjjs.end() ) pdfDirs.erase( predIt++ );
    else ++predIt;
  }
  // Set up PDF to run
  Prediction * pdf = new Prediction( pdfDirs, data );

  Statitical_Effect * statEff = new Statitical_Effect( pdf->pdfFit( "PredictionFunctionForError" ),
						       pdf->covarianceMaticies() );

  //---------------------------------------------------------------------------
  SetAtlasStyle();
  foreach( const double& s, scales ) plotErrors( data, *pdf , *statEff, s );

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


void plotErrors( const Experiment& exp, const Prediction& pdf, const Statitical_Effect& err, const double& l ) {

  TFile outF( (g_outDir+"errors.root").c_str(), "UPDATE" );

  TCanvas c( "c", "", 500, 500 );
  c.cd();
  c.SetLogx();
  gStyle->SetOptStat(0);

  TH1D dummy("dummy","",2,0.9,33);
  dummy.SetXTitle("#chi");

  const int nobs = 11;
  const double xbins[nobs+1] = { 1.0, 1.34986, 1.82212, 2.4596, 3.32012, 4.48169, 6.04965, 8.16617, 11.0232, 14.8797,
                                 20.0855, 30.0 };

  // for each mjj bin
  foreach( const double& mjj, exp.mjjs() )
  {
    string name = str( boost::format( "Err_m%2.0f_L%2.1f" ) % mjj % l );
    string title = str( boost::format( "Errors m_{jj} > %2.0f GeV at #Lambda = %2.1f" ) % mjj % l );

    TH1D errorHisto( name.c_str(), title.c_str(), nobs, xbins );
    errorHisto.SetFillColor(1);
    errorHisto.SetFillStyle(3003);
    
    //pdf.setMjj( mjj );
    // for each chi bin
    for( int bin = 0; bin < exp[mjj].chi().size(); ++bin ) {

      double chi = pdf.labelChi( exp[mjj].chi( bin ) );

      // Evaluate function at scale
      double n = pdf( mjj, chi,  pow( l, -4 ) ); // modify by systematics

      // evaluate errors at scale
      double e = err.error( n, l, chi, mjj );

      errorHisto.SetBinContent( bin+1, n );
      errorHisto.SetBinError( bin+1, e );
    }

    errorHisto.Scale( 1. , "width" );

    dummy.SetMinimum( 0.8 * errorHisto.GetMinimum() );
    dummy.SetMaximum( 1.2 * errorHisto.GetMaximum() );

    dummy.Draw();
    errorHisto.Draw( "SAMEE2" );
    c.Print( (g_outDir+name+".pdf").c_str() );
    c.Print( (g_outDir+name+".png").c_str() );

    name = "rel_"+name;
    title = "Relative "+title;
    TH1D relErrHisto( name.c_str(), title.c_str(), nobs, xbins );

    for ( int bin = 1; bin <= errorHisto.GetNbinsX(); ++bin ) {
      double y = errorHisto.GetBinContent( bin );
      double e = errorHisto.GetBinError( bin );
      relErrHisto.SetBinContent( bin, e / y );
    }

    title = str( boost::format( "d#sigma/#sigma [Full MC #Lambda = %2.1f TeV]" ) % l );
    dummy.SetYTitle( title.c_str() );
    dummy.SetMinimum( 0. );
    dummy.SetMaximum( 1.5 * relErrHisto.GetMaximum() );
    dummy.Draw();
    relErrHisto.Draw( "SAME][" );
    c.Print( (g_outDir+name+".pdf").c_str() );
    c.Print( (g_outDir+name+".png").c_str() );

    outF.cd();
    errorHisto.Write();
    relErrHisto.Write();
    
  }

  
  outF.Close();

  return;

}
