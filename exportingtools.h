#ifndef EXPORTINGTOOLS_H
#define EXPORTINGTOOLS_H

#include <AIS_InteractiveContext.hxx>
#include <meshdatabase.h>
#include <TopoDS_Compound.hxx>

class exportingTools
{
public:

    exportingTools();
    static void exportCloud(const occHandle(AIS_InteractiveContext) &theCTX, meshDataBase *theDB);
    static void exportSTL(const occHandle(AIS_InteractiveContext) &theCTX, meshDataBase *theDB);
    static bool exportSTEP(const TopoDS_Compound &aComp);
};

#endif // EXPORTINGTOOLS_H
