#ifndef SMOOTHINGTOOLS_H
#define SMOOTHINGTOOLS_H

//! ---
//! Qt
//! ---
#include <QMap>

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include <ng_meshvs_datasourceface.h>

class smoothingTools
{
public:

    //! smooth a scalar field
    static void scalarFieldSmoother(QMap<int,double> &field,
                                    const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                    const std::map<int,double> &betaAveField,
                                    const std::map<int,int> &mapOfNodeTypes,
                                    int NbSteps);

    //! smooth a vectorial field
    static void fieldSmoother(QMap<int,QList<double>> &field,
                              const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                              const std::map<int,double> &betaAveField,
                              const std::map<int,double> &betaVisibilityField,
                              const std::map<int,int> &mapOfNodeTypes,
                              int NbSteps,
                              bool normalize);

    //! helper
    static void surroundingNodes(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, int nodeID, bool isLocal, std::set<int> &setOfSurroundingNodes);
};

#endif // SMOOTHINGTOOLS_H
