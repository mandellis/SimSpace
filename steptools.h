#ifndef STEPTOOLS_H
#define STEPTOOLS_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

#include <STEPControl_Reader.hxx>
#include <TColStd_ListOfAsciiString.hxx>

class STEPTools
{
public:

    STEPTools(){;}
    static void dumpSTEPLabels(const STEPControl_Reader &reader, TColStd_ListOfAsciiString &listOfNames);
};

#endif // STEPTOOLS_H
