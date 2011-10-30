#include "TestStatMonitor.hpp"
#include <cmath>
#include <string>

#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "TH1.h"
#include "TH2.h"

using namespace std;
using namespace boost::assign;

TestStatMonitor::TestStatMonitor()
  : _folder("figures/")
  , _ext   (".png") {
  init();
}


TestStatMonitor::TestStatMonitor( const string& folder, const string& ext )
  : _folder( folder )
  , _ext   ( ext ) {
  init();
}


void
TestStatMonitor::init() {

  vector<double> alphas = list_of(0)             (1/pow(8000,2)) (1/pow(7000,2)) (1/pow(6000,2))
                                 (1/pow(5000,2)) (1/pow(4000,2)) (1/pow(3000,2)) (1/pow(1500,2))
                                 (1/pow(1000,2)) (1/pow( 750,2)) (1/pow( 500,2));
  foreach( double alpha, alphas ) {
    _likelihoods     [alpha] = new TH2D("","",100,0,-1,100,0,-1);
    _likelihoodRatios[alpha] = new TH2D("","",100,0,-1,100,0,-1);
    _minimizedAlpha  [alpha] = new TH1D("","",100,alphas.front(),alphas.back());
  }
}


TestStatMonitor::~TestStatMonitor() {
}


void
TestStatMonitor::monitor( Likelihood& l ) {

  

}


void
TestStatMonitor::monitor( LikelihoodRatio& launda ) {
}

