#include "PseudoExperiment.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <stdexcept>

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/checked_delete.hpp>
#include <boost/lexical_cast.hpp>

#include <TCanvas.h>
#include <TGraph.h>
#include <TH2.h>

#include "Prediction.hpp"

using namespace std;
using namespace boost;

struct adder: public unary_function< double, void > {
  adder() :
      sum( 0 ) {
  }
  double sum;
  void operator()( double x ) {
    sum += x;
  }
};

string RandomString( const int len = 10 ) {
  static string charset = "0123456789"
                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                          "abcdefghijklmnopqrstuvwxyz";
  string result( len, 'a' );
  for( int i = 0; i < len; i++ )
    result[i] = charset[rand() % charset.size()];
  return result;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

Experiment::Experiment() {

}

//_____________________________________________________________________________________________________________________
Experiment::Experiment( const map< double, MjjExperiment > subExperiments ) :
    _subExperiments( subExperiments ) {

  typedef const map< double, MjjExperiment > expMap_t;
  foreach( const expMap_t::value_type& exp, _subExperiments )
    _mjjs.insert( exp.first );

}

//_____________________________________________________________________________________________________________________
Experiment::Experiment( const map< double, TH1* > histos ) {
  typedef const map< double, MjjExperiment > expMap_t;
  typedef const map< double, TH1* > histMap_t;
  foreach( const histMap_t::value_type& hist, histos )
  {
    expMap_t::value_type entry = make_pair( hist.first, MjjExperiment( *( hist.second ) ) );
    _subExperiments.insert( entry );
    _mjjs.insert( hist.first );
  }
}

//_____________________________________________________________________________________________________________________
Experiment::~Experiment() {
}

//_____________________________________________________________________________________________________________________
ostream& operator<<( ostream& out, const Experiment& e ) {
  typedef const map< double, MjjExperiment > expMap_t;
  foreach( const expMap_t::value_type& exp, e._subExperiments )
    out << " mjj = " << exp.first << "\n" << exp.second;
  return out;
}

//_____________________________________________________________________________________________________________________
const MjjExperiment& Experiment::mjjExperiment( const double& mjj ) const {
  map< double, MjjExperiment >::const_iterator itr = _subExperiments.find( mjj );
  if ( itr == _subExperiments.end() ) throw( range_error(
      "in " + string( __FILE__ ) + " line " + lexical_cast< string >( __LINE__ ) + ": No experiment for mjj = "
      + lexical_cast< string >( mjj ) ) );
  return itr->second;
}

const MjjExperiment& Experiment::operator[]( const double& mjj ) const {
  return mjjExperiment( mjj );
}

//_____________________________________________________________________________________________________________________
set< double > Experiment::mjjs() const {
  return _mjjs;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

MjjExperiment::MjjExperiment() :
    _integral( 0 ),
    _canvas( 0 ),
    _graph( 0 ) {
}

MjjExperiment::MjjExperiment( const TH1& h ) :
    _name( h.GetName() ),
    _canvas( 0 ),
    _graph( 0 ) {

  _integral = h.Integral();
  int nBins = h.GetNbinsX();
  for( int bin = 1; bin <= nBins; ++bin ) {
    _chi.push_back( h.GetBinCenter( bin ) );
    _n.push_back( h.GetBinContent( bin ) );
  }

}

MjjExperiment::MjjExperiment( const vector< double >& x, const vector< double >& y ) :
    _chi( x ),
    _n( y ),
    _canvas( 0 ),
    _graph( 0 ) {
  _integral = for_each( y.begin(), y.end(), adder() ).sum;
}

MjjExperiment::~MjjExperiment() {
  //delete _canvas;
  //delete _graph;
}

string MjjExperiment::name() const {
  return _name;
}
double MjjExperiment::chi( int& i ) const {
  return _chi.at( i );
}
double MjjExperiment::n( int& i ) const {
  return _n.at( i );
}
vector< double > MjjExperiment::chi() const {
  return _chi;
}
vector< double > MjjExperiment::n() const {
  return _n;
}
double MjjExperiment::integral() const {
  return _integral;
}

void MjjExperiment::name( const string& name ) {
  _name = name;
}
void MjjExperiment::chi( const std::vector< double >& x ) {
  _chi = x;
}
void MjjExperiment::n( const std::vector< double >& y ) {
  _n = y;
  _integral = for_each( y.begin(), y.end(), adder() ).sum;
}

void MjjExperiment::plot() const {

  if ( _canvas || _graph ) {
    cout << "You already plotted this experiment ( " + _name + " ) silly!" << endl;
    return;
  }

  string name;
  if ( _name.empty() ) name = RandomString(); // random string
  else name = _name;
  _canvas = new TCanvas( ( name + "Canvas" ).c_str(), "", 500, 500 );
  _canvas->cd();
  _graph = new TGraph( _chi.size(), &_chi[0], &_n[0] );
  _graph->Draw( "AP" );
  _graph->GetXaxis()->SetTitle( "#chi" );
  _graph->GetYaxis()->SetTitle( _name.c_str() );

}

void MjjExperiment::print( ostream& out ) const {

  out << "Experiment " << _name << '\n';
  assert( _chi.size() == _n.size() );
  for( int i = 0; i < _chi.size(); ++i ) {
    out << "   chi = " << _chi.at( i ) << ",   n = " << _n.at( i ) << '\n';
  }
  out << " sum = " << _integral << '\n';

}
ostream& operator<<( ostream& out, const MjjExperiment& e ) {

  e.print( out );
  return out;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

PseudoExperiment::PseudoExperiment() :
    _alpha( -1. ) {
}

PseudoExperiment::PseudoExperiment( const std::map< double, MjjExperiment >& subExperiments, const double& alpha ) :
    Experiment( subExperiments ),
    _alpha( alpha ) {
}

PseudoExperiment::PseudoExperiment( const std::map< double, TH1* >& hists, const double& alpha ) :
    Experiment( hists ),
    _alpha( alpha ) {
}

double PseudoExperiment::alpha() const {
  return _alpha;
}

void PseudoExperiment::print( ostream& out ) const {

  out << "PEs with alpha = " << _alpha << '\n';
  typedef const map< double, MjjExperiment > expMap_t;
  foreach( const expMap_t::value_type& exp, _subExperiments )
    out << "mjj = " << exp.first << "\n" << exp.second << "\n";

}

ostream& operator<<( ostream& out, const PseudoExperiment& e ) {

  e.print( out );
  return out;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

PseudoExperimentFactory::PseudoExperimentFactory() :
    _pdf( 0 ) {
}

PseudoExperimentFactory::PseudoExperimentFactory( const Prediction* pdf, const Experiment& graft, unsigned int seed ) :
    _pdf( pdf ),
    _graft( graft ),
    _random( seed ) {
}

PseudoExperimentFactory::~PseudoExperimentFactory() {
}

PseudoExperiment PseudoExperimentFactory::build( const double& alpha ) {

  ++_nGenerated[alpha];
  _pdf->newPE();

  map< double, MjjExperiment > mjjPEs;
  vector< double > content;
  vector< double > chis;
  foreach( const double& mjj, _graft.mjjs() )
  {
    string peName = str( format( "PE_mjj%2.1f_L%2.1e_n%1.0f" ) % mjj % pow( alpha, -0.25 ) % _nGenerated[alpha] );
    for( int bin = 0; bin < _graft[mjj].chi().size(); ++bin ) {
      double chi = _graft[mjj].chi( bin );
      double expectedN = ( *_pdf )( mjj, chi, alpha ); // modify by systematics
      chis.push_back( chi );
      content.push_back( _random.Poisson( expectedN ) );
    }

    pair< map< double, MjjExperiment >::iterator, bool > val = mjjPEs.insert(
        make_pair( mjj, MjjExperiment( chis, content ) ) );
    if ( !val.second ) throw( logic_error(
        string( __FILE__ ) + " line " + lexical_cast< string >( __LINE__ ) + ": Generated two PEs at same m_jj = "
        + lexical_cast< string >( mjj ) ) );
    val.first->second.name( peName );
  }

  PseudoExperiment result( mjjPEs, alpha );
  return result;

}

vector< PseudoExperiment > PseudoExperimentFactory::build( const double& alpha, const int& n ) {

  vector< PseudoExperiment > result;
  result.reserve( n );
  for( int i = 0; i < n; ++i )
    result.push_back( build( alpha ) );

  return result;

}

const Prediction * PseudoExperimentFactory::pdf() const {
  return _pdf;
}

