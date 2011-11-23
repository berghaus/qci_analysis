/*
 * RandomFloat.hpp
 *
 *  Created on: 2011-11-18
 *      Author: frank
 */

#ifndef RANDOMFLOAT_HPP_
#define RANDOMFLOAT_HPP_

#include <cassert>
#include <cstdlib>

template< class T >
class Random {
  T _min;
  T _max;
  T _range;
public:
  Random() :
      _min( 0 ),
      _max( 1 ),
      _range( 1 ) {
  }

  Random( const T& min, const T& max ) :
      _min( min ),
      _max( max ),
      _range( max - min ) {
  }
  virtual ~Random() {
  }
  ;

  T operator()() {
    assert( _min < _max );
    T random = ( (T) rand() ) / (T) RAND_MAX;
    return _min + _range * random;
  }

  T operator()( const T& min, const T& max ) {
    _min = min;
    _max = max;
    _range = max - min;
    return ( *this )();
  }

};

#endif /* RANDOMFLOAT_HPP_ */
