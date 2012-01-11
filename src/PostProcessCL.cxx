#include "PostProcessCL.hpp"
/*
 * PostProcessCL.cxx
 *
 *  Created on: 2012-01-11
 *      Author: frank
 */

using namespace std;

PostProcessCL::PostProcessCL() :
    _doCLs( false ) {
}

PostProcessCL::PostProcessCL( const vector< PValueTest >& sigLL ) :
    _doCLs( false ),
    _sigLL( sigLL ) {
}

PostProcessCL::PostProcessCL( const vector< PValueTest >& sigLL, const vector< PValueTest >& bkgLL ) :
    _doCLs( sigLL.size() == bkgLL.size() ),
    _sigLL( sigLL ),
    _bkgLL( bkgLL ) {

}

PostProcessCL::~PostProcessCL() {

}

