#ifndef EXPORTINGTOOLS_H
#define EXPORTINGTOOLS_H

#include <AIS_InteractiveContext.hxx>
#include <meshdatabase.h>
#include <TopoDS_Compound.hxx>

class SimulationNodeClass;

class exportingTools
{
public:

    exportingTools();
    static void exportCloud(const occHandle(AIS_InteractiveContext) &theCTX, meshDataBase *theDB);
    static void exportSTL(const occHandle(AIS_InteractiveContext) &theCTX, meshDataBase *theDB);
    static bool exportSTEP(const TopoDS_Compound &aComp);

    static void exportNodalResult(SimulationNodeClass *aNode, const string &fileName);
};

#endif // EXPORTINGTOOLS_H
