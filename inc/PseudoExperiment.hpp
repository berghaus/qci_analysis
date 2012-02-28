#ifndef PSEUDO_EXPERIMENT_HPP
#define PSEUDO_EXPERIMENT_HPP
#include <set>
#include <ostream>
#include <map>
#include <vector>

#include <TH1.h>
#include <TRandom3.h>

class TCanvas;
class TGraph;

class MjjExperiment;
class Prediction;

class Experiment {
public:
  Experiment();
  Experiment( const std::map< double, MjjExperiment > );
  Experiment( const std::map< double, TH1* > );
  virtual ~Experiment();
  virtual const MjjExperiment& mjjExperiment( const double& ) const;
  virtual std::set< double > mjjs() const;
  virtual const MjjExperiment& operator[] ( const double& ) const;
  friend std::ostream& operator<<( std::ostream&, const Experiment& );
protected:
  std::map< double, MjjExperiment > _subExperiments;
  std::set< double > _mjjs;

};

class MjjExperiment {

public:
  MjjExperiment();
  MjjExperiment( const TH1& );
  MjjExperiment( const std::vector< double >&, const std::vector< double >& );
  virtual ~MjjExperiment();

  double chi( int& ) const;
  double n( int& ) const;

  std::string name() const;
  std::vector< double > chi() const;
  std::vector< double > n() const;
  double integral() const;

  void name( const std::string& );
  void chi( const std::vector< double >& );
  void n( const std::vector< double >& );

  virtual void plot() const;

  friend std::ostream& operator<<( std::ostream&, const MjjExperiment& );

protected:
  virtual void print( std::ostream& ) const;

private:
  std::string _name;
  std::vector< double > _chi;
  std::vector< double > _n;
  double _integral;
  mutable TCanvas *_canvas;
  mutable TGraph *_graph;

};

class PseudoExperiment: public Experiment {

public:
  PseudoExperiment();
  PseudoExperiment( const std::map< double, MjjExperiment >&, const double& );
  PseudoExperiment( const std::map< double, TH1* >&, const double& );
  double alpha() const;

  friend std::ostream& operator<<( std::ostream&, const PseudoExperiment& );

protected:
  virtual void print( std::ostream& ) const;

private:
  double _alpha;

};

class PseudoExperimentFactory {

public:
  PseudoExperimentFactory( const Prediction* pdf, const Experiment& graft, unsigned int seed = 65539 );
  virtual ~PseudoExperimentFactory();
  std::vector< PseudoExperiment > build( const double& alpha, const int& n );
  PseudoExperiment build( const double& alpha );

  // getters
  const Prediction * pdf() const;

private:
  PseudoExperimentFactory();

  const Prediction * _pdf;
  const Experiment _graft;

  TRandom3 _random;
  std::map< double, unsigned int > _nGenerated;
  std::map< double, std::vector< PseudoExperiment* > > _generated;

};

#endif // PSEUDO_EXPERIMENT_HPP
