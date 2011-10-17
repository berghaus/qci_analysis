#ifndef PDF_HPP
#define PDF_HPP

#include <vector>

class TH2;

class PDF {

public:
  PDF();
  PDF( const TH2* );
  PDF( const PDF& );
  virtual ~PDF();

  double operator() ( const double&, const int&, const std::vector<double>& ) const;

  TH2* hist() const;
  void hist( const TH2* );

private:
  TH2* _hist;

};


#endif // PDF_HPP
