#include "PseudoExperiment.hpp"
#define BOOST_TEST_MODULE PseudoExperimentTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <TH1.h>

BOOST_AUTO_TEST_CASE( constructors_test )
{
  Experiment exp0;
  BOOST_CHECK_EQUAL( exp0.integral(), int(0) );
  BOOST_CHECK_EQUAL( exp0.name(), std::string("") );
  BOOST_CHECK( exp0.x().empty() );
  BOOST_CHECK( exp0.y().empty() );


  std::vector<double> x;
  for( int i = 0; i < 10; ++i ) x.push_back( i );
  std::vector<double> y( x.size(), 1 );
  Experiment exp1( x, y );
  BOOST_CHECK_EQUAL( exp1.integral(), double( x.size() ) );
  BOOST_CHECK_EQUAL( exp1.name(), std::string("") );
  BOOST_CHECK_EQUAL( exp1.x().size(), x.size() );
  BOOST_CHECK_EQUAL( exp1.y().size(), y.size() );

  TH1D histogram("testHistogram","",x.size()-1, &x[0] );
  for( int i = 0; i < x.size() -1; ++i ) histogram.Fill( x[i], 1. );
  Experiment exp2( histogram );
  BOOST_CHECK_EQUAL( exp2.integral(), double( histogram.Integral() ) );
  BOOST_CHECK_EQUAL( exp2.name(), std::string("") );
  BOOST_CHECK_EQUAL( exp2.x().size(), histogram.GetNbinsX() );
  BOOST_CHECK_EQUAL( exp2.y().size(), histogram.GetNbinsX() );

}

BOOST_AUTO_TEST_CASE( setter_and_getter_test ){
  Experiment exp0;
  std::vector<double> x;
  for( int i = 0; i < 10; ++i ) x.push_back( i );
  std::vector<double> y( x.size(), 1 );

  exp0.x( x );
  BOOST_REQUIRE_EQUAL( exp0.x().size(), x.size() );
  for( int i = 0; i < x.size() ; ++i )
    BOOST_CHECK_EQUAL( exp0.x()[i], x[i] );

  exp0.y( y );
  BOOST_CHECK_EQUAL( exp0.integral(), double( y.size() ) );
  BOOST_REQUIRE_EQUAL( exp0.y().size(), y.size() );
  for( int i = 0; i < y.size() ; ++i )
    BOOST_CHECK_EQUAL( exp0.y()[i], y[i] );
}

BOOST_AUTO_TEST_CASE( plot_test ){
  std::vector<double> x;
  for( int i = 0; i < 10; ++i ) x.push_back( i );
  std::vector<double> y( x.size(), 1 );
  Experiment exp0( x, y );
  exp0.plot();
}

// EOF
