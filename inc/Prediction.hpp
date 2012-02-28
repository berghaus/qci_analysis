#ifndef PREDICTION_HPP
#define PREDICTION_HPP

#include <map>
#include <string>
#include <vector>

#include <TMatrixTSym.h>
#include <TRandom3.h>

#include "Effect.hpp"
#include "PseudoExperiment.hpp"

class TDirectoryFile;
class TF1;
class TGraphErrors;
class PredictionMonitor;

//! Implements MC prediction of di-jet spectrum in chi for QCI and QCD
/*! \class Prediction
 *
 * Stores the result of MC analysis on QCD and QCI samples with various contact
 * interaction scales (Lambda).  For each bin in chi the prediction is fit with
 *
 *   mu(Lambda) = [0] + [1] * x + [2] * sqrt( x )
 *
 *     where x = Lambda^-4
 *
 * The prediction in each chi bin is then based of these fits to provide well
 * behaved behaviour in Lambda.  The fit is motivated by the contact
 * interaction term added to the standard model Lagrangian.
 */
class Prediction {

public:

  //! default constructor - only use if you know what you are doing
  Prediction();

  //! main constructor - should only be called once TGraph::Fit leaks memory:
  /*!
   * FIXME: use some variation of singelton here to ensure single call behaviour
   * \param a file organised containing one graph per data bin with the
   *        predicted number of events as a function of the scale parameter
   * \param number of events in data to scale the MC prediction to
   */
  Prediction( std::map< double, TDirectoryFile* >, const Experiment& );

  //! recommended constructor - use this once you know the fit parameters
  /*!
   * \param fit parameters
   * \param number of events in data to scale the MC prediction to
   */
  Prediction( const std::map< double, std::vector< double > >&, double );

  //! copy constructor - probably easier to use than recommended constructor
  /*!
   * When using this constructor remember to set the number of events in your
   * experiment after
   * \param Prediction to copy, should have been initialised with main constructor
   */
  Prediction( const Prediction& );

  //! the destructor
  virtual ~Prediction();

  //! for use by the log likelihood class
  /*!
   * scales prediction to nData and compares given number of events to
   * Prediction at chi and CI scale
   * \param di-jet mass (m_jj)
   * \param di-jet separation chi = exp( | y1 - y2 | )
   * \param number of events observed at given chi
   * \param vector with single entry that is the contact interaction parameter
   *        ( Lambda^-4 ) for prediction
   * \return the natural log of the poisson probability
   */
  double operator()( const double&, const double&, const int&, const std::vector< double >& ) const;

  //! expected number of events at chi and CI parameter
  /*!
   * provides predicted number of events not scaled to nData for use in
   * creation of pseudo-experiments
   * \param di-jet mass (m_jj)
   * \param di-jet separation chi = exp( | y1 - y2 | )
   * \param contact interaction parameter ( Lambda^-4 ) for prediction
   * \return expected number of events at chi and alpha modified by
   *         statistical and systematic uncertainty
   */
  double operator()( const double&, const double&, const double& ) const;

  //! currently does nothing
  void init();

  //! accept monitoring visitor
  void accept( PredictionMonitor& );

  //! provides predicted number of events scaled to nData
  /*!
   * scaling divides sum of predicted events at each chi bin at the given CI
   * scale and multiplies by nData
   * \param di-jet separation chi = exp( | y1 - y2 | )
   * \param contact interaction parameter ( Lambda^-4 )
   * \return number of events for given chi and CI scaled to data
   */
  double interpolate( const double&, const double& ) const;

  //! getter from graphs read in by Prediction
  std::map< double, TGraphErrors* > eventCounts( const double& ) const;

  //! getter for parameters resulting from fitting to eventCounts
  std::map< double, std::vector< double > > pdfFitParams( const double& ) const;

  //! return a clone of the pdf fit function
  /*!
   * clones the fit function and returns the pointer. passes ownership of
   * cloned function
   * \param name for TF1 clone
   */
  TF1* pdfFit( const std::string& ) const;

  //! getter for covariance matricies
  std::map< double, TMatrixTSym< double > > covarianceMaticies( const double& ) const;

  //! getter for covariance matricies
  std::map< double, std::map< double, TMatrixTSym< double > > > covarianceMaticies() const;

  //! sum event prediction over chi bins in the given CI parameter
  /*!
   * evaluates fit at given CI scale parameter for each chi and returns result
   * \param contact interaction parameter ( Lambda^-4 ) for prediction
   * \return sum of event prediction in each chi bin
   */
  double sumOverChi( const double& ) const;

  //! setter for number of observed events to normalize the prediction to
  void nData( const Experiment& );

  //! getter for number of observed events to normalize the prediction to
  int nData( const double& ) const;

  //! setter function for statistical and systematics effects to include
  void effects( std::vector< Effect* > );

  //! add a single effect
  void addEffect( Effect* );

  //! set up systematic effects for a new PE
  void newPE() const;

  friend class PredictionMonitor;

private:

  //! fit function to use - same for all chi bins
  TF1* _pdfFit;

  std::set<double> _mjjs;

  class MjjPrediction {
  public:
    MjjPrediction() {
    }
    MjjPrediction( TDirectoryFile*, const double, TF1* );
    virtual ~MjjPrediction() {
    }

    //! number of observed events to normalize the prediction to
    double _nData;

    //! graphs read in containing number of predicted events vs CI parameter
    std::map< double, TGraphErrors* > _eventCounts;

    //! parameters of pdfFit resulting from fitting to event counts
    std::map< double, std::vector< double > > _pdfFitParams;

    //! covariance of pdfFit results
    std::map< double, TMatrixTSym< double > > _covarianceMaticies;

  private:
    TF1* _pdfFit;

  };

  // necessary parameters for each mjj
  std::map< double, MjjPrediction > _predictions;

  //! graphs read in containing number of predicted events vs CI parameter
  //std::map< double, std::map< double, TGraphErrors* > > _eventCounts;

  //! parameters of pdfFit resulting from fitting to event counts
  //std::map< double, std::map< double, std::vector< double > > > _pdfFitParams;

  //! covariance of pdfFit results
  //std::map< double, std::map< double, TMatrixTSym< double > > > _covarianceMaticies;

  //! translate given chi to one used as label in the maps
  double labelChi( const double& ) const;

  //! vector of statistical and systematic errors
  std::vector< Effect* > _effects; // TODO: make smart pointers

  void setMjj( const double& ) const;

  //! mjj for current bin
  mutable double _mjj;

  //! prediction for current bin
  mutable MjjPrediction * _mjjPred;

};


#endif // PREDICTION_HPP
