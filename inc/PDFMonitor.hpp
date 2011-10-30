#ifndef PDF_MONITOR_HPP
#define PDF_MONITOR_HPP

#include <string>

class PDF;

class PDFMonitor {

public:
  PDFMonitor();
  PDFMonitor( const std::string&, const std::string );
  ~PDFMonitor();

  virtual void monitor( PDF& );

private:

  std::string _folder;
  std::string _ext;

};


#endif // PDF_MONITOR_HPP
