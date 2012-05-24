#include "PseudoExperiment.hpp"
#define BOOST_TEST_MODULE PseudoExperimentTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include "PDF.hpp"

BOOST_AUTO_TEST_CASE( constructors_test ) {
  MjjExperiment exp0;
  BOOST_CHECK_EQUAL( exp0.integral(), int(0) );
  BOOST_CHECK_EQUAL( exp0.name(), std::string("") );
  BOOST_CHECK( exp0.chi().empty() );
  BOOST_CHECK( exp0.n().empty() );

  std::vector< double > x;
  for( int i = 0; i < 10; ++i )
    x.push_back( i );
  std::vector< double > y( x.size(), 1 );
  MjjExperiment exp1( x, y );
  BOOST_CHECK_EQUAL( exp1.integral(), double( x.size() ) );
  BOOST_CHECK_EQUAL( exp1.name(), std::string("") );
  BOOST_CHECK_EQUAL( exp1.chi().size(), x.size() );
  BOOST_CHECK_EQUAL( exp1.n().size(), y.size() );

  TH1D histogram( "testHistogram", "", x.size() - 1, &x[0] );
  for( int i = 0; i < x.size() - 1; ++i )
    histogram.Fill( x[i], 1. );
  MjjExperiment exp2( histogram );
  BOOST_CHECK_EQUAL( exp2.integral(), double( histogram.Integral() ) );
  BOOST_CHECK_EQUAL( exp2.name(), std::string(histogram.GetName()) );
  BOOST_CHECK_EQUAL( exp2.chi().size(), histogram.GetNbinsX() );
  BOOST_CHECK_EQUAL( exp2.n().size(), histogram.GetNbinsX() );

  PseudoExperiment pe0;
  BOOST_CHECK_EQUAL( pe0.integral(), int(0) );
  BOOST_CHECK_EQUAL( pe0.name(), std::string("") );
  BOOST_CHECK( pe0.chi().empty() );
  BOOST_CHECK( pe0.n().empty() );
  BOOST_CHECK_EQUAL( pe0.alpha(), -1. );

  double a = 5.;
  PseudoExperiment pe1( x, y, a );
  BOOST_CHECK_EQUAL( pe1.integral(), double( x.size() ) );
  BOOST_CHECK_EQUAL( pe1.name(), std::string("") );
  BOOST_CHECK_EQUAL( pe1.chi().size(), x.size() );
  BOOST_CHECK_EQUAL( pe1.n().size(), y.size() );
  BOOST_CHECK_EQUAL( pe1.alpha(), a );

  TFile * dataF = TFile::Open( "~/docs/comp/analysis/data.root", "READ" );
  TH1 * dataHist = (TH1*) dataF->Get( "Chi_2000-to-7000all" );
  MjjExperiment atlas( *dataHist );
  TFile * input = TFile::Open( "~/docs/comp/analysis/vanilla.root", "READ" );
  Prediction pdf( input, atlas.integral() );
  // int nPE = 10;
  // PseudoExperimentFactory pef0( &pdf, atlas );

  input->Close();
  dataF->Close();

}

BOOST_AUTO_TEST_CASE( setter_and_getter_test ) {
  MjjExperiment exp0;
  std::vector< double > x;
  for( int i = 0; i < 10; ++i )
    x.push_back( i );
  std::vector< double > y( x.size(), 1 );

  exp0.chi( x );
  BOOST_REQUIRE_EQUAL( exp0.chi().size(), x.size() );
  for( int i = 0; i < x.size(); ++i )
    BOOST_CHECK_EQUAL( exp0.chi()[i], x[i] );

  exp0.n( y );
  BOOST_CHECK_EQUAL( exp0.integral(), double( y.size() ) );
  BOOST_REQUIRE_EQUAL( exp0.n().size(), y.size() );
  for( int i = 0; i < y.size(); ++i )
    BOOST_CHECK_EQUAL( exp0.n()[i], y[i] );

  PseudoExperiment pe0;
  pe0.chi( x );
  BOOST_REQUIRE_EQUAL( pe0.chi().size(), x.size() );
  for( int i = 0; i < x.size(); ++i )
    BOOST_CHECK_EQUAL( pe0.chi()[i], x[i] );

  pe0.n( y );
  BOOST_CHECK_EQUAL( pe0.integral(), double( y.size() ) );
  BOOST_REQUIRE_EQUAL( pe0.n().size(), y.size() );
  for( int i = 0; i < y.size(); ++i )
    BOOST_CHECK_EQUAL( pe0.n()[i], y[i] );

}

BOOST_AUTO_TEST_CASE( plot_test ) {
  std::vector< double > x;
  for( int i = 0; i < 10; ++i )
    x.push_back( i );
  std::vector< double > y( x.size(), 1 );
  MjjExperiment exp0( x, y );
  exp0.plot();

  double a = 5.;
  PseudoExperiment pe0( x, y, a );
  pe0.plot();

}

// EOF
