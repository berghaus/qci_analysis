/*
 * Error.cxx
 *
 *  Created on: 2012-02-12
 *      Author: frank
 */

#include "Error.hpp"

#include <TF1.h>

using namespace std;

Statitical_Error::Statitical_Error() :
    _fitFunction( 0 ) {
}

Statitical_Error::Statitical_Error( const TF1* fitFunction,
                                    const map< double, TMatrixTSym< double > >& covarianceMaticies ) :
    _fitFunction( (TF1*) fitFunction->Clone( "fitFunction" ) ),
    _covarianceMaticies( covarianceMaticies ) {

}

Statitical_Error::~Statitical_Error() {
  delete _fitFunction;
}

double Statitical_Error::apply( const double& mu, const double& lambda, const double& chi, const double& mjj ) {

  double result = mu;

  // do stuff to result.

  return result;
}

