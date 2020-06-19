#ifndef SURFACEMESHCUTTER_H
#define SURFACEMESHCUTTER_H

#include "occhandle.h"
#include <ng_meshvs_datasourceface.h>

class surfaceMeshCutter
{
public:

    surfaceMeshCutter();

    static bool cutSurfaceMesh(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax,
                               const occHandle(MeshVS_DataSource) &inputMeshDS,
                               occHandle(MeshVS_DataSource) &cutMesh);

    static bool cutSurfaceMesh(const occHandle(MeshVS_DataSource) &inputMeshDS,
                               const QList<TopoDS_Shape> &listOfShapes,
                               occHandle(MeshVS_DataSource) &cutMesh);

};

#endif // SURFACEMESHCUTTER_H
