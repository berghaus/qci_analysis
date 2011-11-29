#ifndef PSEUDO_EXPERIMENT_HPP
#define PSEUDO_EXPERIMENT_HPP
#include <map>
#include <vector>

#include <TH1.h>
#include <TRandom.h>

#include "PDF.hpp"

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

  void plot();

private:
  std::string _name;
  std::vector<double> _x;
  std::vector<double> _y;
  double _integral;

};


class PseudoExperiment : public Experiment {

public:
  PseudoExperiment( const std::vector<double>&, const std::vector<double>&, const double& );
  double alpha() const;

private:
  std::vector<double> _x;
  std::vector<double> _y;
  double _integral;
  double _alpha;

};


class PseudoExperimentFactory {

public:
  PseudoExperimentFactory( const PDF* pdf, const Experiment& graft, unsigned int seed = 65539 );
  virtual ~PseudoExperimentFactory();
  std::vector<PseudoExperiment> build( const double& alpha, const int& n );

private:

  PseudoExperimentFactory();
  const PDF *      _pdf;
  const Experiment _graft;

  TRandom _random;
  std::map<double,unsigned int>      _nGenerated;
  std::map<double,std::vector<PseudoExperiment*> > _generated;

  PseudoExperiment build( const double& alpha );

};


#endif // PSEUDO_EXPERIMENT_HPP
