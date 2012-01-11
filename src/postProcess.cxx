/*
 * postProcess.cxx
 *
 *  Created on: 2012-01-11
 *      Author: frank
 */
#include <string>

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

#define foreach BOOST_FOREACH
#define ERROR_NO_SIGNAL_INPUT 1

using namespace std;
namespace po = boost::program_options;

int main( int argc, char* argv[] ) {

  // process cmd opts
  // Declare the supported options.
  po::options_description desc( "Allowed options" );
  desc.add_options()( "help,h", "print this help message" )
      ( "sigInputFiles,s", po::value< vector< string > >()->multitoken(), "list of input files containing signal likelihood distributions" )
      ( "bkgInputFiles,b", po::value< vector< string > >()->multitoken(), "list of input files containing background likelihood distributions" )
      ( "outDir,o", po::value< string >(), "output directory for plots" );

  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, desc ), vm );
  po::notify( vm );

  if ( vm.count( "help" ) ) {
    cout << desc << "\n";
    return 1;
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

  string folder = "./";
  if ( vm.count( "outDir" ) ) {
    folder = vm["outDir"].as<string>();
  }
  cout << "directing figures to: " << folder << "\n";

  return 0;
}

