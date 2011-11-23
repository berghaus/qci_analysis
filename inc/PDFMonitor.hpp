#ifndef PDF_MONITOR_HPP
#define PDF_MONITOR_HPP

#include <string>
#include <vector>
#include <TCanvas.h>


class PDF;
class TGraph;

class PDFMonitor {

public:
  PDFMonitor();
  PDFMonitor( const std::string&, const std::string );
  virtual ~PDFMonitor();

  virtual void monitor( PDF& );

private:

  std::string _folder;
  std::string _ext;
  std::vector<TGraph*> _interpolations;
  std::vector<TGraph*> _fitResults;

  TCanvas _pdfCanvas;
  TCanvas _interpolCanvas;
  TCanvas _fitResultCanvas;


};


#endif // PDF_MONITOR_HPP
