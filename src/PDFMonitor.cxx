#include "PDFMonitor.hpp"

#include "PDF.hpp"
using namespace std;


PDFMonitor::PDFMonitor()
  : _folder ( "figures/PDF/" )
  , _ext    ( ".png" ) {
}


PDFMonitor::PDFMonitor( const string& folder, const string ext )
  : _folder( folder )
  , _ext   ( ext ) {
}


PDFMonitor::~PDFMonitor() {
}


void
PDFMonitor::monitor( PDF& pdf ) {

  // make some histograms
  pdf;

}

