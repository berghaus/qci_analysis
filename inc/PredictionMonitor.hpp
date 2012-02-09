#ifndef PREDICTION_MONITOR_HPP
#define PREDICTION_MONITOR_HPP

#include <string>
#include <vector>
#include <TCanvas.h>

class Prediction;
class TGraph;

class PredictionMonitor {

public:
  PredictionMonitor();
  PredictionMonitor( const std::string&, const std::string );
  virtual ~PredictionMonitor();

  virtual void monitor( Prediction& );

private:

  std::string _folder;
  std::string _ext;
  std::vector<TGraph*> _interpolations;
  std::vector<TGraph*> _fitResults;
  std::vector<TGraph*> _parameters;

  TCanvas _pdfCanvas;
  TCanvas _interpolCanvas;
  TCanvas _fitResultCanvas;
  TCanvas _parameterCanvas;


};


#endif // PREDICTION_MONITOR_HPP
