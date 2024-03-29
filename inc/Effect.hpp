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
class TH1;
class TFile;
namespace appl {
class grid;
}
;

//! Interface for systematic and stochastic errors
/*! \class Error
 * interface class for stochastic and systematic effects to follow
 */
class Effect {
public:
  virtual ~Effect() {
  }

  virtual double apply( double, const double&, const double&, const double& ) const = 0;
  virtual void newPE() = 0; //< select new _nSigma

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
  Statitical_Effect( const TF1*, const std::map< double, std::map< double, TMatrixTSym< double > > >& );
  virtual ~Statitical_Effect();

  //! implementation of function from Effect interface
  virtual double apply( double, const double&, const double&, const double& ) const;
  virtual void newPE() {
  } //< nothing needs to be done here
  double error( const double&, const double&, const double&, const double& ) const;

private:

  TF1 * _fitFunction;
  std::map< double, std::map< double, TMatrixTSym< double > > > _covarianceMaticies;

  mutable TRandom3 _random;

};

//! Implements systematic error descriped as function of chi and mjj
class Experimental_Systematic_Effect: public Effect {
public:
  Experimental_Systematic_Effect();
  Experimental_Systematic_Effect( const std::string& );
  virtual ~Experimental_Systematic_Effect();

  //! implementation of function from Effect interface
  virtual double apply( double, const double&, const double&, const double& ) const;

protected:
  double _nSigma;
private:
  TFile * _file;
  std::map< double, std::map< double, std::map< double, TH1* > > > _errors;
};

class JES_Effect: public Experimental_Systematic_Effect {
public:
  JES_Effect();
  JES_Effect( const std::string& );
  virtual ~JES_Effect();
  virtual void newPE(); //< select new _nSigma
private:
  TRandom3 _random;
};

class JER_Effect: public Experimental_Systematic_Effect {
public:
  JER_Effect();
  JER_Effect( const std::string& );
  virtual ~JER_Effect();
  virtual void newPE(); //< select new _nSigma
private:
  TRandom3 _random;
};

class PDF_Effect: public Effect {
public:
  PDF_Effect();
  PDF_Effect( const std::string& );
  virtual ~PDF_Effect();

  virtual double apply( double, const double&, const double&, const double& ) const;
  virtual void newPE(); //< select new _nSigma
private:
  double _nSigma;
  TRandom3 _random;
  TFile * _file;
  std::map< double, std::pair< double, double > > _errors; //< key = label chi, val = + / - error
};


class Scale_Effect: public Effect {
public:
  Scale_Effect();
  Scale_Effect( const std::vector< std::string >& );
  virtual ~Scale_Effect();

  virtual double apply( double, const double&, const double&, const double& ) const;
  virtual void newPE(); //< select new _nSigma
private:
  int _nLoops;
  double _eScale;
  double _muf;
  double _mur;
  TRandom3 _random;
  appl::grid * _grid;
  TH1 * _nominal;
  TH1 * _changed;
  TF1 * _distribution;
};

#endif // EFFECT_HPP
