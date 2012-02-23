#include "PseudoExperiment.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/checked_delete.hpp>

#include <TCanvas.h>
#include <TGraph.h>
#include <TH2.h>

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

Experiment::Experiment() :
    _integral( 0 ),
    _canvas( 0 ),
    _graph( 0 ) {
}

Experiment::Experiment( const TH1& h ) :
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

Experiment::Experiment( const vector< double >& x, const vector< double >& y ) :
    _chi( x ),
    _n( y ),
    _canvas( 0 ),
    _graph( 0 ) {
  _integral = for_each( y.begin(), y.end(), adder() ).sum;
}

Experiment::~Experiment() {
  //delete _canvas;
  //delete _graph;
}

string Experiment::name() const {
  return _name;
}
double Experiment::chi( int& i ) const {
  return _chi.at( i );
}
double Experiment::n( int& i ) const {
  return _n.at( i );
}
vector< double > Experiment::chi() const {
  return _chi;
}
vector< double > Experiment::n() const {
  return _n;
}
double Experiment::integral() const {
  return _integral;
}

void Experiment::name( const string& name ) {
  _name = name;
}
void Experiment::chi( const std::vector< double >& x ) {
  _chi = x;
}
void Experiment::n( const std::vector< double >& y ) {
  _n = y;
  _integral = for_each( y.begin(), y.end(), adder() ).sum;
}

void Experiment::plot() const {

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

void Experiment::print( ostream& out ) const {

  out << "Experiment " << _name << '\n';
  assert( _chi.size() == _n.size() );
  for( int i = 0; i < _chi.size(); ++i ) {
    out << "   chi = " << _chi.at( i ) << ",   n = " << _n.at( i ) << '\n';
  }
  out << " sum = " << _integral << '\n';

}
ostream& operator<<( ostream& out, const Experiment& e ) {

  e.print( out );
  return out;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

PseudoExperiment::PseudoExperiment() :
    _alpha( -1. ) {
}

PseudoExperiment::PseudoExperiment( const vector< double >& x, const vector< double >& y, const double& alpha ) :
    Experiment( x, y ),
    _alpha( alpha ) {
}

double PseudoExperiment::alpha() const {
  return _alpha;
}

void PseudoExperiment::print( ostream& out ) const {

  out << "Experiment " << name() << '\n';
  assert( chi().size() == n().size() );
  for( int i = 0; i < chi().size(); ++i ) {
    out << "   chi = " << chi( i ) << ",   n = " << n( i ) << '\n';
  }
  out << " sum =   " << integral() << '\n';
  out << " alpha = " << _alpha << '\n';

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

  string peName = str( format( "PE_alpha%2.1e_n%1.0f" ) % pow( alpha, -0.25 ) % _nGenerated[alpha] );

  vector< double > content;
  vector< double > chis;
  for( int bin = 0; bin < _graft.chi().size(); ++bin ) {
    double chi = _graft.chi( bin );
    double expectedN = ( *_pdf )( chi, alpha ); // modify by systematics
    chis.push_back( chi );
    content.push_back( _random.Poisson( expectedN ) );
  }

  PseudoExperiment result( chis, content, alpha );
  result.name( peName );

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

