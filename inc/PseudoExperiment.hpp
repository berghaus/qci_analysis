#ifndef PSEUDO_EXPERIMENT_HPP
#define PSEUDO_EXPERIMENT_HPP
#include <ostream>
#include <map>
#include <vector>

#include <TH1.h>
#include <TRandom.h>

#include "PDF.hpp"

class TCanvas;
class TGraph;



class Experiment {

public:
  Experiment();
  Experiment( const TH1& );
  Experiment( const std::vector<double>&, const std::vector<double>& );
  virtual ~Experiment();

  double x( int& ) const;
  double y( int& ) const;

  std::string name() const;
  std::vector<double> x() const;
  std::vector<double> y() const;
  double integral() const;

  void name( const std::string& );
  void x( const std::vector<double>& );
  void y( const std::vector<double>& );

  virtual void plot() const;

  friend std::ostream& operator<< ( std::ostream&, const Experiment& );

protected:
  virtual void print( std::ostream& ) const;

private:
  std::string _name;
  std::vector<double> _x;
  std::vector<double> _y;
  double _integral;
  mutable TCanvas *_canvas;
  mutable TGraph *_graph;

};


class PseudoExperiment : public Experiment {

public:
  PseudoExperiment();
  PseudoExperiment( const std::vector<double>&, const std::vector<double>&, const double& );
  double alpha() const;

  friend std::ostream& operator<< ( std::ostream&, const PseudoExperiment& );

protected:
  virtual void print( std::ostream& ) const;

private:
  double _alpha;

};


class PseudoExperimentFactory {

public:
  PseudoExperimentFactory( const PDF* pdf, const Experiment& graft, unsigned int seed = 65539 );
  virtual ~PseudoExperimentFactory();
  std::vector<PseudoExperiment> build( const double& alpha, const int& n );
  PseudoExperiment build( const double& alpha );

  // getters
  const PDF * pdf() const;

private:
  PseudoExperimentFactory();

  const PDF *      _pdf;
  const Experiment _graft;

  TRandom _random;
  std::map<double,unsigned int>      _nGenerated;
  std::map<double,std::vector<PseudoExperiment*> > _generated;

};


#endif // PSEUDO_EXPERIMENT_HPP
