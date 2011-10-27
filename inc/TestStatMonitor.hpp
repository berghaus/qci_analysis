#ifndef TEST_STAT_MONITOR_HPP
#define TEST_STAT_MONITOR_HPP
#include <map>
#include <string>

class TH1;
class TH2;
class Likelihood;
class LikelihoodRatio;

class TestStatMonitor {

public:

  TestStatMonitor();
  TestStatMonitor( const std::string&, const std::string& );
  virtual ~TestStatMonitor();

  virtual void monitor( Likelihood& );
  virtual void monitor( LikelihoodRatio& );

private:
  void init();

  std::string _folder;
  std::string _ext;

  std::map<double,TH2*> _likelihoods;
  std::map<double,TH2*> _likelihoodRatios;
  std::map<double,TH1*> _minimizedAlpha;

};


#endif // TEST_STAT_MONITOR_HPP
