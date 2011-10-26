#ifndef TEST_STAT_MONITOR_HPP
#define TEST_STAT_MONITOR_HPP
#include <string>


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
  std::string _folder;
  std::string _ext;

};


#endif // TEST_STAT_MONITOR_HPP
