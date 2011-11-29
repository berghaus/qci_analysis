#include "PseudoExperiment.hpp"
#include <algorithm>
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


struct adder : public unary_function<double, void> {
  adder() : sum(0) {}
  double sum;
  void operator()(double x) { sum += x; }
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

Experiment::Experiment() {
}


Experiment::Experiment( const TH1& h ) {

  _integral = h.Integral();
  int nBins = h.GetNbinsX();
  for( int bin = 1; bin <= nBins; ++bin ) {
    _x.push_back( h.GetBinCenter( bin ) );
    _y.push_back( h.GetBinContent( bin ) );
  }

}


Experiment::Experiment( const vector<double>& x, const vector<double>& y ) :
  _x( x ),
  _y( y ) {
  _integral = for_each( y.begin(), y.end(), adder() ).sum;
}


Experiment::~Experiment() {
}

string Experiment::name() const { return _name; }
double Experiment::x( int& i ) const { return _x.at(i); }
double Experiment::y( int& i ) const { return _y.at(i); }
vector<double> Experiment::x() const { return _x; }
vector<double> Experiment::y() const { return _y; }
double Experiment::integral() const { return _integral; }

void Experiment::name( const string& name ) { _name = name; }
void Experiment::x( const std::vector<double>& x) {  _x = x; }
void Experiment::y( const std::vector<double>& y) {
  _y = y;
  _integral = for_each( y.begin(), y.end(), adder() ).sum;
}

void Experiment::plot() {

  TCanvas * dataC = new TCanvas( (_name+"Canvas").c_str(), "", 500, 500 );
  dataC->cd();
  TGraph * g = new TGraph( _x.size(), &_x[0], &_y[0] );
  g->Draw("AP");
  g->GetXaxis()->SetTitle( "#chi" );
  g->GetYaxis()->SetTitle( _name.c_str() );

}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

PseudoExperiment::PseudoExperiment( const vector<double>& x, const vector<double>& y, const double& alpha ) :
  Experiment( x, y ),
  _alpha( alpha ) {
  }

double PseudoExperiment::alpha() const { return _alpha; }



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



PseudoExperimentFactory::PseudoExperimentFactory( const PDF* pdf, const Experiment& graft, unsigned int seed ) :
  _pdf( pdf ),
  _graft( graft ),
  _random( seed ) {
}

PseudoExperimentFactory::~PseudoExperimentFactory() {
}

PseudoExperiment
PseudoExperimentFactory::build( const double& alpha ) {

  ++_nGenerated[alpha];

  string peName = str( format( "PE_alpha%2.1e_n%1.0f" ) % pow( alpha, -0.25 ) % _nGenerated[alpha] );

  vector<double> content;
  vector<double> chis;
  for( int bin = 0; bin < _graft.x().size(); ++bin ) {
    double chi = _graft.x( bin );
    double expectedN = ( *_pdf )( chi, alpha );
    chis.push_back( chi );
    content.push_back( _random.Poisson( expectedN ) );
  }

  PseudoExperiment result( chis, content, alpha );
  result.name( peName );

  // cout << "Integral( " << peName << " ) = " << result.integral() << '\n'; 
  return result;

}


vector< PseudoExperiment > PseudoExperimentFactory::build( const double& alpha, const int& n ) {

  vector< PseudoExperiment > result;
  result.reserve( n );
  for( int i = 0; i < n; ++i ) result.push_back( build( alpha ) );

  return result;

}


PseudoExperimentFactory::PseudoExperimentFactory() :
  _pdf( 0 ) {
}
