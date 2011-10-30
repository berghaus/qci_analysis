#include "PDFMonitor.hpp"
#include <vector>
#include <boost/format.hpp>


#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"


#include "PDF.hpp"
using namespace std;
using boost::format;


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
  TCanvas c("PDFMonCanvas","",600,800); c.cd();

  for( int xBin = 1; xBin <= pdf.hist()->GetNbinsX(); ++xBin ) {
    string projName = string( pdf.hist()->GetName() ) + str( format("-xbin-%1.0i" ) % xBin );
    TH1* proj = pdf.hist()->ProjectionY( projName.c_str(), xBin, xBin );

    proj->Draw("EX0");
    c.Print( (_folder+projName+_ext).c_str() );
  }

}

