#ifndef PDF_HPP
#define PDF_HPP

#include <map>
#include <vector>

class TFile;
class TGraphErrors;
class PDFMonitor;

class PDF {

public:
  PDF();
  PDF( TFile* );
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

  std::map<double,TGraphErrors*> eventCounts() const;

private:
  TFile* _file;
  std::map<double,TGraphErrors*> _eventCounts;

};


#endif // PDF_HPP
