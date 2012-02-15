/*
 * Error.cxx
 *
 *  Created on: 2012-02-12
 *      Author: frank
 */

#include "Error.hpp"

#include <cmath>
#include <stdexcept>

#include <boost/assign.hpp>
#include <boost/format.hpp>

#include <TF1.h>

using boost::format;

using namespace std;
using namespace boost::assign;

//_____________________________________________________________________________________________________________________
Statitical_Error::Statitical_Error() :
    _fitFunction( 0 ) {
}

//_____________________________________________________________________________________________________________________
Statitical_Error::Statitical_Error( const TF1* fitFunction,
                                    const map< double, TMatrixTSym< double > >& covarianceMaticies ) :
    _fitFunction( (TF1*) fitFunction->Clone( "fitFunction" ) ),
    _covarianceMaticies( covarianceMaticies ) {

}

//_____________________________________________________________________________________________________________________
Statitical_Error::~Statitical_Error() {
  delete _fitFunction;
}

//_____________________________________________________________________________________________________________________
double Statitical_Error::apply( const double& mu, const double& lambda, const double& chi, const double& mjj ) {

  double result = mu;
  double err = error( mu, lambda, chi, mjj );

  // predicted number of events modified by statistical uncertainty
  result = _random.Gaus( result, err ); // must be positive

  return result > 0. ?
      result :
      0.;

}

//_____________________________________________________________________________________________________________________
double Statitical_Error::error( const double& mu, const double& lambda, const double& chi, const double& mjj ) {

  double alpha = pow( lambda, -4. );

  // do stuff to result.
  if ( _covarianceMaticies.find( chi ) == _covarianceMaticies.end() ) 
    throw( range_error( str( format( "Statistical_Error::error - no covariance matrix at chi = %2.1f in " + string( __FILE__ ) + ": %3.0i" ) % chi % __LINE__ ) ) );

  const TMatrixTSym< double > & mat = _covarianceMaticies.find( chi )->second;

  // I know the gradients of my fit function are:
  vector< double > grad = list_of( 0. )( alpha )( pow( alpha, 0.5 ) );
  vector< double > sum( 3, 0. );
  double error = 0.;

  // first sum
  for( int i = 0; i < mat.GetNrows(); ++i ) {
    sum[i] = 0;
    for( int j = 0; j < mat.GetNcols(); ++j ) {
      sum[i] += mat( i, j ) * grad[j];
    }
  }

  // second sum
  for( int i = 0; i < mat.GetNcols(); ++i ) {
    error += sum[i] * grad[i];
  }

  return error;

}
