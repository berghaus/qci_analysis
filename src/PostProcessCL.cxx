#include "PostProcessCL.hpp"
/*
 * PostProcessCL.cxx
 *
 *  Created on: 2012-01-11
 *      Author: frank
 */
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "CertaintyLevel.hpp"
#include "Likelihood.hpp"

using namespace std;

PostProcessCL::PostProcessCL() :
    _doCLs( false ),
    _CLsb( 0 ),
    _CLs( 0 ) {
}

PostProcessCL::PostProcessCL( const vector< PValueTest >& sigLL, vector< Neg2LogLikelihoodRatio* >& errorLLRs, Neg2LogLikelihoodRatio* dataLLR ) :
    _doCLs( false ),
    _sigLL( sigLL ),
    _errorLLRs( errorLLRs ),
    _dataLLR( dataLLR ),
    _CLsb( 0 ),
    _CLs( 0 ) {
}

PostProcessCL::PostProcessCL( const vector< PValueTest >& sigLL, const vector< PValueTest >& bkgLL, vector< Neg2LogLikelihoodRatio* >& errorLLRs, Neg2LogLikelihoodRatio* dataLLR ) :
    _doCLs( sigLL.size() == bkgLL.size() ),
    _sigLL( sigLL ),
    _bkgLL( bkgLL ),
    _errorLLRs( errorLLRs ),
    _dataLLR( dataLLR ),
    _CLsb( 0 ),
    _CLs( 0 ) {

}

PostProcessCL::~PostProcessCL() {

  if( _CLsb ) delete _CLsb;
  if( _CLs ) delete _CLs;

}

// ensure PValue vectors are sorted in alpha
bool compByAlpha( const PValueTest& x, const PValueTest& y ) {
  return x.alpha() > y.alpha(); // want largest alpha (= lowest lambda) first
}

void PostProcessCL::proc() {

  sort( _sigLL.begin(), _sigLL.end(), compByAlpha );
  if ( _doCLs ) sort( _bkgLL.begin(), _bkgLL.end(), compByAlpha );

  _CLsb = new CertaintyLevel( "CL_{s+b}", _sigLL.size(), pow( _sigLL.front().alpha(), -0.25 ),
                              pow( _sigLL.back().alpha(), -0.25 ) );
  if ( _doCLs ) _CLs = new CertaintyLevel( "CL_{s}", _sigLL.size(), pow( _sigLL.front().alpha(), -0.25 ),
                                           pow( _sigLL.back().alpha(), -0.25 ) );

  vector< PValueTest >::iterator sigItr = _sigLL.begin();
  vector< PValueTest >::iterator sigEnd = _sigLL.end();
  vector< PValueTest >::iterator bkgItr = _bkgLL.begin();
  vector< PValueTest >::iterator bkgEnd = _bkgLL.end();

  while ( sigItr != sigEnd && ( !_doCLs || bkgItr != bkgEnd ) ) {

    if ( _doCLs && sigItr->alpha() != bkgItr->alpha() ) {
      cout << "error mismatch in singal and background likelihood distributions!\n";
      continue;
    }

    double alpha = sigItr->alpha();
    double scale = pow( alpha, -0.25 );

    vector< double > par( 1, alpha );
    double clsb_observed = (*sigItr)( *_dataLLR );

    vector< double > clsb_expected;
    clsb_expected.reserve( _errorLLRs.size() );
    foreach( Neg2LogLikelihoodRatio * l, _errorLLRs )
    clsb_expected.push_back( (*sigItr)( *l ) );
    sort( clsb_expected.begin(), clsb_expected.end() );

    _CLsb->add( scale, clsb_observed, clsb_expected );

    if ( _doCLs ) {

      double cls_observed = (*sigItr)( *_dataLLR ) / (*bkgItr)( *_dataLLR );

      vector< double > cls_expected;
      cls_expected.reserve( _errorLLRs.size() );
      foreach( Neg2LogLikelihoodRatio* l, _errorLLRs )
        cls_expected.push_back( (*sigItr)( *l ) / (*bkgItr)( *l ) );
      sort( cls_expected.begin(), cls_expected.end() );

      _CLs->add( scale, cls_observed, cls_expected );
    }

    ++sigItr;
    if ( _doCLs ) ++bkgItr;
  }

}


void PostProcessCL::print() {

  if ( _CLsb ) cout << *_CLsb << "\n\n";
  if ( _CLs ) cout << *_CLs << "\n\n";

}


void PostProcessCL::plot( const string& folder ) {

  if ( _CLsb ) _CLsb->plot( folder );
  if ( _CLs ) _CLs->plot( folder );

}


