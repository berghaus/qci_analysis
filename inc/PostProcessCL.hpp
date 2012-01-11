/*
 * PostProcessCL.hpp
 *
 *  Created on: 2012-01-11
 *      Author: frank
 */

#ifndef POSTPROCESSCL_HPP_
#define POSTPROCESSCL_HPP_

#include <vector>

#include "PValueTest.hpp"

class PostProcessCL {

public:
  PostProcessCL();
  PostProcessCL( const std::vector< PValueTest >& );
  PostProcessCL( const std::vector< PValueTest >&, const std::vector< PValueTest >& );
  virtual ~PostProcessCL();

private:

  bool _doCLs;

  std::vector< PValueTest > _sigLL;
  std::vector< PValueTest > _bkgLL;

};


#endif /* POSTPROCESSCL_HPP_ */
