/*
 * Error.hpp
 *
 *  Created on: 2012-02-12
 *      Author: frank
 */
#ifndef ERROR_HPP
#define ERROR_HPP
#include <map>

#include <TMatrixTSym.h>
#include <TRandom3.h>

class TF1;

// error application interface
struct Error {

  Error() {
  }
  virtual ~Error() {
  }

  virtual double apply( const double&, const double&, const double&, const double& ) = 0;

};

// Error implementation for statical effect
class Statitical_Error: public Error {
public:
  Statitical_Error();
  Statitical_Error( const TF1*, const std::map< double, TMatrixTSym< double > >& );
  virtual ~Statitical_Error();

  virtual double apply( const double&, const double&, const double&, const double& );

private:

  double error( const double&, const double&, const double&, const double& );
  TF1 * _fitFunction;
  std::map< double, TMatrixTSym< double > > _covarianceMaticies;

  TRandom3 _random;

};

#endif // ERROR_HPP
