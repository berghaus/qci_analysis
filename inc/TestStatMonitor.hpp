#ifndef TEST_STAT_MONITOR_HPP
#define TEST_STAT_MONITOR_HPP
#include <map>
#include <string>

class TH1;
class TH2;
class Likelihood_FCN;
class LikelihoodRatio_FCN;

class TestStatMonitor {

public:

  TestStatMonitor();
  TestStatMonitor( const std::string&, const std::string& );
  virtual ~TestStatMonitor();

  virtual void monitor( Likelihood_FCN& );
  virtual void monitor( LikelihoodRatio_FCN& );

  virtual void finalize();

private:
  void init();

  std::string _folder;
  std::string _ext;

  TH2* _likelihood;
  TH2* _likelihoodRatio;
  TH1* _minimizedAlpha;
  TH1* _minimizedLaunda;

};

#endif // TEST_STAT_MONITOR_HPP
