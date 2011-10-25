#ifndef PDF_HPP
#define PDF_HPP

#include <vector>

#include <TProfile2D.h>

class TGraph2D;
class TH2;

class PDF {

public:
  PDF();
  PDF( TH2* );
  PDF( const PDF& );
  virtual ~PDF();

  // returns the natural log of the poisson probability of number of events
  // (2nd argument) using the provided x (1st parameter) and parameters (3rd)
  // to look up the number of events epxected by the PDF.  Interpolates for
  // smoothness
  double operator() ( const double&, const int&, const std::vector<double>& ) const;

  // expectated number of events at x, y
  double operator() ( const double&, const double& ) const;

  TH2* hist() const;
  void hist( const TH2* );

private:
  TH2* _hist;
  TGraph2D *_graph;

  TProfile2D* _monitor;

};


#endif // PDF_HPP
