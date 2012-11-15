/*
 * PostProcessCL.hpp
 *
 *  Created on: 2012-01-11
 *      Author: frank
 */

#ifndef POSTPROCESSCL_HPP_
#define POSTPROCESSCL_HPP_

#include <vector>
#include <string>

#include "PValueTest.hpp"

class Neg2LogSimpleLikelihoodRatio;
class CertaintyLevel;

class PostProcessCL {

public:
  PostProcessCL();
  // Constructor for CLs+b limit only.  Does not rake ownership of likelihood ratios supplied
  PostProcessCL( const std::vector< PValueTest<Neg2LogSimpleLikelihoodRatio> >&, std::vector< Neg2LogSimpleLikelihoodRatio* >&, Neg2LogSimpleLikelihoodRatio* );

  // Constructor for CLs+b and CLs limit.  Does not rake ownership of likelihood ratios supplied
  PostProcessCL( const std::vector< PValueTest<Neg2LogSimpleLikelihoodRatio> >&, const std::vector< PValueTest<Neg2LogSimpleLikelihoodRatio> >&,
                 std::vector< Neg2LogSimpleLikelihoodRatio* >&, Neg2LogSimpleLikelihoodRatio* );
  virtual ~PostProcessCL();

  virtual void proc();
  virtual void print();
  virtual void plot( const std::string& folder = "./" );

private:

  bool _doCLs;

  std::vector< PValueTest<Neg2LogSimpleLikelihoodRatio> > _sigLL;
  std::vector< PValueTest<Neg2LogSimpleLikelihoodRatio> > _bkgLL;
  std::vector< Neg2LogSimpleLikelihoodRatio* > _errorLLRs;
  Neg2LogSimpleLikelihoodRatio* _dataLLR;

  double _pValue;
  CertaintyLevel * _CLsb;
  CertaintyLevel * _CLs;

};

#endif /* POSTPROCESSCL_HPP_ */
