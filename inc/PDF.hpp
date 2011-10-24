#ifndef PDF_HPP
#define PDF_HPP

#include <vector>

class TGraph2D;
class TH2;

class PDF {

public:
  PDF();
  PDF( TH2* );
  PDF( const PDF& );
  virtual ~PDF();

   // Probability of observind number of events at x given parameters
  double operator() ( const double&, const int&, const std::vector<double>& ) const;
   // expectated number of events at x, y
  double operator() ( const double&, const double& ) const;

  TH2* hist() const;
  void hist( const TH2* );

private:
  TH2* _hist;
  TGraph2D *_graph;

};


#endif // PDF_HPP
