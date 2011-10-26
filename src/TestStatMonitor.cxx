#include "TestStatMonitor.hpp"
#include <string>

using namespace std;

TestStatMonitor::TestStatMonitor()
  : _folder("figures/")
  , _ext   (".png") {
}


TestStatMonitor::TestStatMonitor( const string& folder, const string& ext )
  : _folder( folder )
  , _ext   ( ext ) {
}


TestStatMonitor::~TestStatMonitor() {
}


void
TestStatMonitor::monitor( Likelihood& l ) {
}


void
TestStatMonitor::monitor( LikelihoodRatio& launda ) {
}

