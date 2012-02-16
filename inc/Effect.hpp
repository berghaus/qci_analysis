/*
 * Error.hpp
 *
 *  Created on: 2012-02-12
 *      Author: frank
 */
#ifndef EFFECT_HPP
#define EFFECT_HPP
#include <map>

#include <TMatrixTSym.h>
#include <TRandom3.h>

class TF1;

//! Interface for systematic and stochastic errors
/*! \class Error
 * interface class for stochastic and systematic effects to follow
 */
class Effect {
public:
  virtual ~Effect() {
  }

  virtual double apply( const double&, const double&, const double&, const double& ) const = 0;

};

//! Implements the statistical error on the QCD and QCI MC simulation
/*! \class Statitical_Effect
 * Multiplies the covariance matrix for given chi with the gradient of the
 * Prediction fit function in parameter space to determine the stochastic error
 * on the prediction.  For use when generating pseudo-experiments.
 *
 * Assumes the fit function supplied is
 *
 *     mu(Lambda) = [0] + [1] * x + [2] * sqrt( x )
 *
 *         where x = Lambda^-4
 *
 * Such that the gradient wrt. to the parameters are
 *
 *          dmu/d[0] = 0.
 *          dmu/d[1] = alpha
 *          dmu/d[2] = sqrt( x ) = Lambda^-2
 *
 */
class Statitical_Effect: public Effect {
public:
  Statitical_Effect();
  Statitical_Effect( const TF1*, const std::map< double, TMatrixTSym< double > >& );
  virtual ~Statitical_Effect();

  //! implementation of function from Effect interface
  virtual double apply( const double&, const double&, const double&, const double& ) const;

private:

  double error( const double&, const double&, const double&, const double& ) const;
  TF1 * _fitFunction;
  std::map< double, TMatrixTSym< double > > _covarianceMaticies;

  mutable TRandom3 _random;

};

#endif // EFFECT_HPP
