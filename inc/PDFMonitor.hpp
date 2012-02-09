#ifndef PDF_MONITOR_HPP
#define PDF_MONITOR_HPP

#include <string>
#include <vector>
#include <TCanvas.h>


class Prediction;
class TGraph;

class PDFMonitor {

public:
  PDFMonitor();
  PDFMonitor( const std::string&, const std::string );
  virtual ~PDFMonitor();

  virtual void monitor( Prediction& );

private:

  std::string _folder;
  std::string _ext;
  std::vector<TGraph*> _interpolations;
  std::vector<TGraph*> _fitResults;
  std::vector<TGraph*> _parameters;

  TCanvas _pdfCanvas;
  TCanvas _interpolCanvas;
  TCanvas _fitResultCanvas;
  TCanvas _parameterCanvas;


};


#endif // PDF_MONITOR_HPP
