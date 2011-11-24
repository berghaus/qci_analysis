#ifndef TEST_STAT_MONITOR_HPP
#define TEST_STAT_MONITOR_HPP
#include <map>
#include <string>

#include "Random.hpp"

class TH1;
class TH2;
class Likelihood_FCN;
class LikelihoodRatio;

class TestStatMonitor {

public:

  TestStatMonitor();
  TestStatMonitor( const std::string&, const std::string& );
  virtual ~TestStatMonitor();

  virtual void monitor( Likelihood_FCN& );
  virtual void monitor( LikelihoodRatio& );

  virtual void finalize();

private:
  void init();

  std::string _folder;
  std::string _ext;

  TH2* _likelihoodVsScale;
  TH2* _likelihoodRatioVsScale;
  TH2* _likelihoodVsAlpha;
  TH2* _likelihoodRatioVsAlpha;
  TH1* _minimizedAlpha;
  TH1* _minimizedLaunda;

  Random<double> _randomCompScale;
  Random<double> _randomAlpha;

  double _minScale;
  double _maxScale;
  int _nBinsScale;

};

#endif // TEST_STAT_MONITOR_HPP
