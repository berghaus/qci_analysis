#ifndef CERTAINTY_LEVEL_HPP
#define CERTAINTY_LEVEL_HPP
#include <ostream>
#include <string>
#include <vector>
#include <boost/math/distributions/normal.hpp>

template< class T >
T quantile( const std::vector< T >&, const double& );

double findX( const double&, const std::vector< double >&, const std::vector< double >& );

class CertaintyLevel {

public:
  // default constructor - should not be used
  CertaintyLevel();

  // constructor specifying name (CLs, CLs+b, etc) and number of values to expect
  CertaintyLevel( const std::string&, const std::vector< double >::size_type& );

  // constructor specifying name (CLs, CLs+b, etc) and number of values to expect and min/max scale to plot
  CertaintyLevel( const std::string&, const std::vector< double >::size_type&, const double&, const double& );

  // adds point in scale with observed and expected with error bands
  void add( const double&, const double&, const std::vector< double >& );

  // draws a canvas with the current graph
  void plot( const std::string& folder = "./" );

  friend std::ostream& operator<<( std::ostream&, const CertaintyLevel& );

private:

  // determine observed 95% exclusion
  double observed() const;
  // determine expected 95% exclusion
  double expected() const;

  int _nBins; // number of scale bins
  double _min; // range for plotting
  double _max; // range for plotting

  std::string _name;

  // TODO: Encapsulate these in a class
  std::vector< double > _scales;
  std::vector< double > _obs;
  std::vector< double > _exp;
  std::vector< double > _e02;
  std::vector< double > _e98;
  std::vector< double > _e16;
  std::vector< double > _e84;
  std::vector< double > _ex;

  boost::math::normal _myGaus;

};

// -----
template< class T >
T quantile( const std::vector< T >& v, const double& q ) {
  typename std::vector< T >::size_type index = v.size() * q;
  return v.at( index );
}

#endif // CERTAINTY_LEVEL_HPP
