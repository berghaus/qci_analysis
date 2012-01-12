#include "PostProcessCL.hpp"
/*
 * PostProcessCL.cxx
 *
 *  Created on: 2012-01-11
 *      Author: frank
 */

#include "Likelihood.hpp"

using namespace std;

PostProcessCL::PostProcessCL() :
    _doCLs( false ) {
}

PostProcessCL::PostProcessCL( const vector< PValueTest >& sigLL, vector< Neg2LogLikelihoodRatio* >& errorLLRs ) :
    _doCLs( false ),
    _sigLL( sigLL ),
    _errorLLRs( errorLLRs ) {
}

PostProcessCL::PostProcessCL( const vector< PValueTest >& sigLL, const vector< PValueTest >& bkgLL, vector< Neg2LogLikelihoodRatio* >& errorLLRs ) :
    _doCLs( sigLL.size() == bkgLL.size() ),
    _sigLL( sigLL ),
    _bkgLL( bkgLL ),
    _errorLLRs( errorLLRs ) {

}

PostProcessCL::~PostProcessCL() {

}

