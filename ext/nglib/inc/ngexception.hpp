#ifndef FILE_NGEXCEPTION
#define FILE_NGEXCEPTION

/**************************************************************************/
/* File:   ngexception.hpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. Jan. 2002                                                  */
/**************************************************************************/
#include <string>
namespace netgen
{

/// Base class for all ng exceptions
class NgException 
{
  /// verbal description of exception
  std::string what;
public:
  ///
  DLL_HEADER NgException (const std::string & s);
  ///
  DLL_HEADER virtual ~NgException ();

  /// append string to description
  //DLL_HEADER void Append (const std::string &s);
    void Append (const char * s);
  
  /// verbal description of exception
  const std::string & What() const { return what; }
};
}

#endif
