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
                                    int NbSteps,
                                    int inflationStep);

    //! smooth a vectorial field
    static void fieldSmoother(QMap<int,QList<double>> &field,
                              const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                              const std::map<int,double> &betaAveField,
                              const std::map<int,double> &betaVisibilityField,
                              const std::map<int,int> &mapOfNodeTypes,
                              int NbSteps,
                              bool normalize,
                              int inflationStep);


    /*
    //! smooth a scalar field - 1
    static void scalarFieldSmoother(QMap<int,double> &field,
                                    const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                    const std::map<int,double> &betaAveField,
                                    const std::map<int,int> &mapOfNodeTypes,
                                    int NbSteps,
                                    int inflationStep);

*/


    //! helpers
    static void surroundingNodes(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, int nodeID, bool isLocal, std::set<int> &setOfSurroundingNodes);
    static bool getPointCoordinate(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS, int globalNodeID, double *P);
    static void smoothScalarAtPoint(const occHandle(Ng_MeshVS_DataSourceFace) &aMeshDS,
                                    const QMap<int,double> &field,
                                    int globalNodeID,
                                    double betaAve,
                                    int inflationStep,
                                    double &fieldValue);

};

#endif // SMOOTHINGTOOLS_H
