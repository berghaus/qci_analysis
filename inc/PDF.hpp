#ifndef PDF_HPP
#define PDF_HPP

#include <map>
#include <string>
#include <vector>

#include <TMatrixTSym.h>
#include <TRandom3.h>

class TDirectoryFile;
class TF1;
class TGraphErrors;
class PDFMonitor;

class PDF {

public:

  // default constructor - should not be used
  PDF();

  // main constructor - should only called ony when needed TGraph::Fit leaks memory:
  // 1st parameter: file organised containing one graph per data bin with the predicted
  //                number of events as a function of the scale parameter
  // 2nd parameter: number of events in data to scale the MC prediction to
  PDF( TDirectoryFile*, const double );

  // recommended constructor:
  // 1st parameter: fit parameters
  // 2nd parameter: number of events in data to scale the MC prediction to
  PDF( const std::map< double, std::vector< double > >&, double );

  // copy constructor - since there are some dynamic member
  PDF( const PDF& );
  virtual ~PDF();

  // returns the natural log of the poisson probability of number of events
  // (2nd argument) using the provided x (1st parameter) and parameters (3rd)
  // to look up the number of events epxected by the PDF.  Interpolates for
  // smoothness
  double operator()( const double&, const int&, const std::vector< double >& ) const;

  // expectated number of events at chi and alpha modified by statistical and systematic uncertainty
  double operator()( const double&, const double& ) const;

  void init();

  void accept( PDFMonitor& );

  double interpolate( const double&, const double& ) const;
  double error( const double&, const double& ) const;

  std::map< double, TGraphErrors* > eventCounts() const;
  std::map< double, std::vector< double > > pdfFitParams() const;
  TF1* pdfFit( const std::string& ) const;
  double sumOverChi( const double& ) const;

  void nData( const int& );
  int nData() const;

private:
  TF1* _pdfFit;
  double _nData;
  std::map< double, TGraphErrors* > _eventCounts;
  std::map< double, std::vector< double > > _pdfFitParams;
  std::map< double, TMatrixTSym<double> > _covarianceMaticies;
  mutable TRandom3 _random;

};

#endif // PDF_HPP
