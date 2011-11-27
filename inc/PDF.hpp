#ifndef PDF_HPP
#define PDF_HPP

#include <map>
#include <string>
#include <vector>

class TFile;
class TF1;
class TGraphErrors;
class PDFMonitor;

class PDF {

public:
  PDF();
  PDF( TFile*, const double );
  PDF( const PDF& );
  virtual ~PDF();

  // returns the natural log of the poisson probability of number of events
  // (2nd argument) using the provided x (1st parameter) and parameters (3rd)
  // to look up the number of events epxected by the PDF.  Interpolates for
  // smoothness
  double operator() ( const double&, const int&, const std::vector<double>& ) const;

  // expectated number of events at x, y
  double operator() ( const double&, const double& ) const;

  void init();

  void accept( PDFMonitor& );

  double interpolate( const double&, const double& ) const;

  std::map< double, TGraphErrors* > eventCounts() const;
  std::map< double, std::vector< double > > pdfFitParams() const;
  TF1* pdfFit( const std::string& ) const;

private:
  TFile* _file;
  TF1* _pdfFit;
  TF1* _normalizedPdfFit;
  double _nData;
  std::map<double,TGraphErrors*> _eventCounts;
  std::map<double,std::vector< double > > _pdfFitParams;

};


#endif // PDF_HPP
