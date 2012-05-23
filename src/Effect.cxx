/*
 * Error.cxx
 *
 *  Created on: 2012-02-12
 *      Author: frank
 */

#include "Effect.hpp"

#include <cmath>
#include <stdexcept>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#ifndef foreach
#define foreach BOOST_FOREACH
#endif
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TKey.h>
#include <TDirectoryFile.h>
#include <TClass.h>

using boost::format;

using namespace std;
using namespace boost::assign;

//=====================================================================================================================
//=====================================    Statitical_Effect   ========================================================
//=====================================================================================================================
//_____________________________________________________________________________________________________________________
Statitical_Effect::Statitical_Effect() :
    _fitFunction( 0 ) {
}

//_____________________________________________________________________________________________________________________
Statitical_Effect::Statitical_Effect( const TF1* fitFunction,
                                      const map< double, map< double, TMatrixTSym< double > > >& covarianceMaticies ) :
    _fitFunction( (TF1*) fitFunction->Clone( "fitFunction" ) ),
    _covarianceMaticies( covarianceMaticies ) {

}

//_____________________________________________________________________________________________________________________
Statitical_Effect::~Statitical_Effect() {
  delete _fitFunction;
}

//_____________________________________________________________________________________________________________________
double Statitical_Effect::apply( const double& mu, const double& lambda, const double& chi, const double& mjj ) const {

//  cout << "in Statitical_Effect::apply\n";
//  cout << "   * mu = " << mu << "\n";
//  cout << "   * lambda = " << lambda << "\n";
//  cout << "   * chi = " << chi << "\n";
//  cout << "   * mjj = " << mjj << "\n";
  double result = mu;
  double err = error( mu, lambda, chi, mjj );
//  cout << "   * err = " << err << "\n";

// predicted number of events modified by statistical uncertainty
  result = _random.Gaus( result, err ); // must be positive

//  cout << "   * modified mu = " << result << "\n";
  return result > 0. ?
      result :
      0.;

}

//_____________________________________________________________________________________________________________________
double Statitical_Effect::error( const double& mu, const double& lambda, const double& chi, const double& mjj ) const {

  double alpha = pow( lambda, -4. );

  // find appropriate covariance matrix
  map< double, map< double, TMatrixTSym< double > > >::const_iterator mjjIt = _covarianceMaticies.lower_bound( mjj );
  if ( mjjIt == _covarianceMaticies.end() ) throw( range_error(
      string( __FILE__ ) +  " "
      + str(
          format( "Statistical_Error::error line %3.0i: - no covariance matrix at mjj = %3.1f in " ) % __LINE__ % mjj ) ) );

  map< double, TMatrixTSym< double > >::const_iterator chiIt = mjjIt->second.find( chi );
  if ( chiIt == mjjIt->second.end() ) throw( range_error(
      string( __FILE__ ) + " "
      + str(
          format( "Statistical_Error::error line %3.0i: - no covariance matrix at chi = %2.1f in " ) % __LINE__ % chi ) ) );

  const TMatrixTSym< double > & mat = chiIt->second;

  // I know the gradients of my fit function are:
  vector< double > grad = list_of( 0. )( alpha )( sqrt( alpha ) );
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



//=====================================================================================================================
//=====================================    Experimental_Systematic_Effect   ========================================================
//=====================================================================================================================

//_____________________________________________________________________________________________________________________
Experimental_Systematic_Effect::Experimental_Systematic_Effect() :
    _file( 0 ),
    _nSigma( 0 ) {
}

//_____________________________________________________________________________________________________________________
Experimental_Systematic_Effect::Experimental_Systematic_Effect( const string& fName ) :
    _file( TFile::Open( fName.c_str(), "READ" ) ),
    _nSigma( 0. ) {

  using boost::lexical_cast;

  // fill our look-up table
  cout << "Reading in systematics from " << fName << "\n";
  // iterate over scale directories
  TIter nextScaleKey( _file->GetListOfKeys() );
  TKey * scaleKey = (TKey*) nextScaleKey();
  do {

    if ( !scaleKey ) break;

    string name = scaleKey->GetName();
    TObject * obj = scaleKey->ReadObj();

    if ( obj && obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      TDirectory * scaleDir = static_cast< TDirectory* >( obj );
      cout << "   * found scale: " << name << endl;
      double scale = lexical_cast< double >( name.substr( 1, name.size() - 1 ) ) / 1000.;
      TIter nextSigmaKey( scaleDir->GetListOfKeys() );
      TKey * sigmaKey = (TKey*) nextSigmaKey();
      do {

        if ( !scaleKey ) break;

        name = sigmaKey->GetName();
        obj = sigmaKey->ReadObj();
        if ( obj && obj->IsA()->InheritsFrom( "TDirectory" ) ) {
          TDirectory * sigmaDir = static_cast< TDirectory* >( obj );
          double sign = name.substr( 0, 1 ) == "p" ?
              1. :
              -1.;
          double sigma = sign * lexical_cast< double >( name.substr( 1, name.size() - 1 ) );
          cout << "      * found sigma: " << sigma << endl;
          TIter nextHistoKey( sigmaDir->GetListOfKeys() );
          TKey * histoKey = (TKey*) nextHistoKey();
          do {

            name = histoKey->GetName();
            obj = histoKey->ReadObj();
            if ( obj && obj->IsA()->InheritsFrom( "TH1" ) ) {
              cout << "         * found histogram: " << name << endl;
              double mjjMin = lexical_cast< double >( name.substr( 0, name.find( "-to-" ) ) );
              double mjjMax = lexical_cast< double >( name.substr( name.find( "-to-" ) + 4, name.size() ) );
              _errors[scale][sigma][mjjMax] = static_cast< TH1* >( obj );
            }

          } while ( histoKey = (TKey*) nextHistoKey() );
        }

      } while ( sigmaKey = (TKey*) nextSigmaKey() );
    }

  } while ( scaleKey = (TKey*) nextScaleKey() );

  // check that all sigma keys are the same distance apart:
  typedef map< double, map< double, map< double, TH1* > > > errorMap;
  double space = 0;
  foreach( errorMap::value_type& sigmas, _errors )
  {
    if ( space == 0 ) {
      space = ( sigmas.second.rbegin()->first - sigmas.second.begin()->first ) / double( sigmas.second.size() - 1 );
    } else if ( space
                != ( sigmas.second.rbegin()->first - sigmas.second.begin()->first ) / double( sigmas.second.size() - 1 ) ) {
      cerr << "WARNING: sigma spacing different at different compositeness scales!\n";
    }
    typedef map< double, map< double, TH1* > > sigmaMap;
    sigmaMap::iterator loItr = sigmas.second.begin();
    sigmaMap::iterator hiItr = sigmas.second.begin();
    ++hiItr;
    sigmaMap::iterator end = sigmas.second.end();
    bool homo = true;
    for( ; hiItr != end && loItr != hiItr; ++loItr, ++hiItr ) {
      if ( space != hiItr->first - loItr->first ) {
        cerr << "WARNING: sigma spacing inhomogenious expect the unexpected!\n";
        homo = false;
        break;
      }
    }
    // if spacing is homogeneous (as it should be) bump everything by half a sapce
    // to let lower_bound find the wanted correction for a given _nSigma
    if ( homo ) {
      sigmaMap::iterator itr = sigmas.second.begin();
      end = sigmas.second.end();
      for( ; itr != end; ++itr ) {
        cout << "changing " << itr->first << " to ";
        const_cast< double& >( itr->first ) = itr->first + space / 2.;
        cout << itr->first << "\n";
      }

    }

  }

}

//_____________________________________________________________________________________________________________________
Experimental_Systematic_Effect::~Experimental_Systematic_Effect() {
}

//_____________________________________________________________________________________________________________________
double Experimental_Systematic_Effect::apply( const double& mu, const double& lambda, const double& chi, const double& mjj ) const {
  double result = mu;

  double lKey = 0.; // <- QCD
  if ( lambda == lambda && lambda < _errors.rbegin()->first ) lKey = _errors.lower_bound( lambda )->first;

  double sKey = 0;
  if ( _nSigma > _errors.find( lKey )->second.rbegin()->first ) sKey = _errors.find( lKey )->second.rbegin()->first;
  else sKey = _errors.find( lKey )->second.lower_bound( _nSigma )->first;

  double mKey = _errors.find( lKey )->second.find( sKey )->second.rbegin()->first;
  if ( mjj < mKey ) mKey = _errors.find( lKey )->second.find( sKey )->second.lower_bound( mjj )->first;

  TH1* histo = _errors.find( lKey )->second.find( sKey )->second.find( mKey )->second;
  double ratio = histo->GetBinContent( histo->FindBin( chi ) );
  result *= ratio;

  return result;
}

//_____________________________________________________________________________________________________________________
JES_Effect::JES_Effect() {

}

//_____________________________________________________________________________________________________________________
JES_Effect::JES_Effect( const std::string& fName ) :
    Experimental_Systematic_Effect( fName ) {

}

//_____________________________________________________________________________________________________________________
JES_Effect::~JES_Effect() {

}

//_____________________________________________________________________________________________________________________
void JES_Effect::newPE() {
  _nSigma = _random.Gaus( 0, 1. );
}

//_____________________________________________________________________________________________________________________
JER_Effect::JER_Effect() {

}

//_____________________________________________________________________________________________________________________
JER_Effect::JER_Effect( const std::string& fName ) :
    Experimental_Systematic_Effect( fName ) {

}

//_____________________________________________________________________________________________________________________
JER_Effect::~JER_Effect() {

}

//_____________________________________________________________________________________________________________________
void JER_Effect::newPE() {
  _nSigma = fabs( _random.Gaus( 0, 1. ) );
}

